// ---------------------------------------------------------
//
//  main.cpp
//  Tyson Brochu 2008
//  Adrian Buganza 2017
//
//  Functions for setting up and running an explicit surface simulation, 
//
// ---------------------------------------------------------


// ---------------------------------------------------------
// Defines
// ---------------------------------------------------------

// Whether to use an OpenGL GUI, or run command-line style
// (Define NO_GUI in your build to suppress the GUI.)

#define NO_GUI

// Whether to run rendering and simulation on separate threads
#define RUN_ASYNC

// Whether to use the El Topo C API.  If not defined, uses the El Topo classes directly.
//#define USE_C_API


// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

// std
#include <cstdio>
#include <fenv_include.h>
#include <fstream>
#include <vector>
#include <queue>

// common
#include <array2.h>
#include <ccd_wrapper.h>
#include <collisionqueries.h>
#include <expansion.h>
#include <marching_tiles_hires.h>
#include <util.h>
#include <vec.h>
#include <wallclocktime.h>

// el topo
#include <collisionpipeline.h>
#include <eltopo.h>
#include <iomesh.h>
#include <meshrenderer.h>
//#include <runstats.h>
#include <surftrack.h>
#include <trianglequality.h>

// UL mechanics solver
#include <framestepper.h>
//#include <meshdriver.h> // only needed for specific examples 
//#include <invagination.h> // only needed for specific examples
#include <scriptinit.h>
#include <simulation.h>

// Eigen, careful
#include <Eigen/Dense>

// diffusion function, skip the include file
void advance_diffusion(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, std::vector<double> &c_A, std::vector<double> &c_B, std::vector<double> &c_C, double dt);

// mechanics (may 17 2018)
void driverVelocities_Explicit(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hfl, std::vector<Vec3d> &A_t,double dt);
void driverVelocities_Uniaxial(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hfl, std::vector<Vec3d> &A_t,double dt, double t_tot);
void driverVelocities_OffBiaxial(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hfl, std::vector<Vec3d> &A_t,double dt, double t_tot);
void driverVelocities_Biaxial(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hfl, std::vector<Vec3d> &A_t,double dt, double t_tot);
void driverVelocities_BulgeTest(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hfl, std::vector<Vec3d> &A_t,double dt, double t_tot);
void updateForces_Explicit(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t, std::vector<Eigen::VectorXd> &ea_t, double dt, double dt_real, const NonDestructiveTriMesh &mesh, double t_tot);
void updateForces_Explicit(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t, std::vector<Eigen::VectorXd> &ea_t, double dt, double dt_real, const NonDestructiveTriMesh &mesh, double t_tot, std::vector<Vec3d> &fa_u, const std::vector<double> &c_C, const std::vector<Vec3d> &x_cell, const NonDestructiveTriMesh &mesh_cell);
void updateForces_Explicit(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t, std::vector<Eigen::VectorXd> &ea_t, 
    std::vector<double> &c_A, std::vector<double> &c_B, std::vector<double> &c_C, std::vector<Vec3d> &fa_u, double dt, double dt_real, const NonDestructiveTriMesh &mesh, double t_tot);

// standalone driver for epibole, january 2019
void driverVelocities_Zebrafish(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt, double t_tot);

#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>    //for windows, you'll need http://www.softagalleria.net/dirent.php or something similar
#include <sys/stat.h>

#ifdef _MSC_VER

#include <pthread.h> //requires use pthreads-win from http://sources.redhat.com/pthreads-win32/, or something like it

//this is kind of hacky, but seems to do the trick for now (on Windows, the name is prefaced by an underscore)
//so we'll just rename is here
#define snprintf _snprintf 

#endif


// ---------------------------------------------------------
// Global interface declarations
// ---------------------------------------------------------

extern "C" {
    void exactinit();    // in predicates.c
}

int main(int argc, char **argv); 

// ---------------------------------------------------------
// Global variable
// ---------------------------------------------------------

//RunStats g_stats;

char g_output_path[256];  // Where to write output data
char g_base_output_path[256];  // Where to write output data
bool g_render_to_texture;

struct hack_camera
{
    float target[3];
    float dist;
    float pitch;
    float heading;
} gluvi_cam;


SurfTrack* g_surf_cell = NULL;        // The dynamic surface for cell
SurfTrack* g_surf_subs = NULL;        // The dynamic surface for the substrate


// Local to main.cpp
namespace {    

    MeshDriver* driver = NULL;     // The thing that makes the dynamic surface go

    // ---------------------------------------------------------
    // Interface declarations for functions local to this file
    // ---------------------------------------------------------
    
    void advance_sim( double dt ); 
    void advance_frame(); 
    void run_simulation( ); 
    void set_up_output_path( );     
    void init_simulation( char *script_filename ); 
    
    // ---------------------------------------------------------
    // Variables local to this file
    // ---------------------------------------------------------
    
    
    Simulation* sim = NULL;        // Timekeeping
    FrameStepper* frame_stepper = NULL;
      
    
    // ---------------------------------------------------------
    ///
    /// Advance simulation, GUI-ignorant
    ///
    // ---------------------------------------------------------
    
    void advance_sim( double sim_dt )
    {
        // sim_dt is the desired time step, but it may not be achieved in a single 
        // step, instead, the step might be split into multiple little steps
        // depending on collisions, self-intersections, weird things happening to 
        // the mesh, etc.
        
        double accum_dt = 0;
        while ( (accum_dt < 0.99 * sim_dt) && (sim->m_curr_t + accum_dt < sim->m_max_t) )
        {            
            std::cout << "\n\n ------- sim step cell ------- \n\n";
            std::cout << "curr t: " << sim->m_curr_t + accum_dt << " / " << sim->m_max_t << std::endl;
            
            std::cout<<"total vertices before improvement = "<<g_surf_cell->get_num_vertices()<<"\n";
            std::cout<<"total triangles before improvemen t= "<<g_surf_cell->m_mesh.num_triangles()<<"\n";
            
            // ---------- 
            // mesh maintenance & topology changes
            // ---------- 
            
            // Improve
            
            g_surf_cell->improve_mesh();
            g_surf_cell->defrag_mesh();
            
            std::cout<<"total vertices after improvement = "<<g_surf_cell->get_num_vertices()<<"\n";
            std::cout<<"total triangles after improvemen t= "<<g_surf_cell->m_mesh.num_triangles()<<"\n";            

            // ---------- 
            // advance underlying simulation
            // ----------
            
            // this is the desired time step for the sub-step, if that makes sense
            double curr_dt = sim_dt - accum_dt; 
            curr_dt = min( curr_dt, sim->m_max_t - sim->m_curr_t - accum_dt );
            
            std::cout << "curr_dt: " << curr_dt << std::endl;

            ////-------------------------------------------------------------////
            //*************** Explicit time integration ***********************//
            ////-------------------------------------------------------------////

            std::vector<Vec3d> new_positions( g_surf_cell->get_num_vertices() );
            // replace the driver with the explicit time integration with UL growth
            //driver->set_predicted_vertex_positions( *g_surf, new_positions, sim->m_curr_t + accum_dt, curr_dt );
            std::vector<Vec3d> V_t_hlf(g_surf_cell->get_num_vertices());
            //std::cout<<"explicit driver...\n";

            // July 13th 2019, off biaxial
            //std::cout<<"explicit driver...\n";
            //driverVelocities_Uniaxial(g_surf->get_positions(), g_surf->get_V_t(), V_t_hlf, g_surf->get_A_t(), curr_dt,sim->m_curr_t );
            //driverVelocities_OffBiaxial(g_surf->get_positions(), g_surf->get_V_t(), V_t_hlf, g_surf->get_A_t(), curr_dt,sim->m_curr_t );
            //driverVelocities_Biaxial(g_surf->get_positions(), g_surf->get_V_t(), V_t_hlf, g_surf->get_A_t(), curr_dt,sim->m_curr_t );
            //driverVelocities_BulgeTest(g_surf->get_positions(), g_surf->get_V_t(), V_t_hlf, g_surf->get_A_t(), curr_dt,sim->m_curr_t );
            driverVelocities_Explicit(g_surf_cell->get_positions(), g_surf_cell->get_V_t(), V_t_hlf, g_surf_cell->get_A_t(), curr_dt);

            // use the V_t_hlf to get desired positions 
            for(int ni=0;ni<g_surf_cell->get_num_vertices();ni++){
                new_positions[ni] = g_surf_cell->get_position(ni)+curr_dt*V_t_hlf[ni];
                //std::cout<<"half velocity "<<V_t_hlf[ni]<<"\n";
                //std::cout<<" new["<<ni<<"]="<<new_positions[ni]<<", from old="<<g_surf->get_position(ni)<<"\n";                
            }
            //std::cout<<"done computing new positions based on V_hlf...\n";
            // this is setting the desired positions
            g_surf_cell->set_all_newpositions( new_positions );
            std::cout << "sim_dt: " << curr_dt << std::endl;
            
            // ----------
            // move & handle collision detection
            // ----------
            
            double actual_dt;
            std::vector<Vec3d> initial_positions = g_surf_cell->get_positions();
            
            // here the dynamic surface would like to indeed take the desired time step
            // but maybe it can't, and in that case it will return the actual time
            // step taken.
            std::cout<<"integrate...\n";
            g_surf_cell->integrate( curr_dt, actual_dt );
            
            // the time step actually taken by el topo
            double ori_dt = curr_dt;
            curr_dt = actual_dt;    
            
            std::cout << "actual_dt: " << actual_dt << std::endl;
            
            std::vector<Vec3d> final_positions = g_surf_cell->get_positions();
            //for(int ni=0;ni<final_positions.size();ni++){
            //    std::cout<<" final["<<ni<<"]="<<final_positions[ni]<<", from initial="<<initial_positions[ni]<<"\n";
            //}
            
            accum_dt += curr_dt;
            
            // with the actual final positions and initial positions recompute the 
            // residual at the real time step, update the velocities, and also the strains
            //updateForces_Explicit(initial_positions, final_positions, g_surf->get_V_t(), V_t_hlf, g_surf->get_A_t(), g_surf->get_ea_t(), ori_dt, actual_dt, g_surf->m_mesh, sim->m_curr_t);
            // update forces including concentrations 
            updateForces_Explicit(initial_positions, final_positions, g_surf_cell->get_V_t(), V_t_hlf, g_surf_cell->get_A_t(), g_surf_cell->get_ea_t(), g_surf_cell->get_cA(),g_surf_cell->get_cB(),g_surf_cell->get_cC(),g_surf_cell->get_fa_u(), ori_dt, actual_dt, g_surf_cell->m_mesh, sim->m_curr_t);
            
            //
            // file output
            //
            
            //char binary_filename[256];
            //sprintf( binary_filename, "%s/frame%04d.bin", g_output_path, frame_stepper->get_frame() );      
            //write_binary_file( g_surf->m_mesh, g_surf->get_positions(), g_surf->m_masses, sim->m_curr_t, binary_filename );   

            //char vtk_filename[256];
            //sprintf( vtk_filename, "%s/frame%04d.vtk", g_output_path, frame_stepper->get_frame() );      
            // several options
            //write_vtk( g_surf->m_mesh, g_surf->get_positions(), sim->m_curr_t, vtk_filename );   
            //write_vtk( g_surf->m_mesh, g_surf->get_positions(), g_surf->get_cB(),g_surf->get_cC(),g_surf->get_cBC(),sim->m_curr_t, vtk_filename );   
            //write_vtk( g_surf->m_mesh, g_surf->get_positions(), g_surf->get_ea_t(), sim->m_curr_t, vtk_filename );   
            
            //char stats_filename[256];
            //sprintf( stats_filename, "%s/aaa-stats.txt", g_output_path );      
            //g_stats.write_to_file( stats_filename );
            
        }
        
        std::cout << " -------------- end sim step cell -------------- \n" << std::endl;


        accum_dt = 0;
        while ( (accum_dt < 0.99 * sim_dt) && (sim->m_curr_t + accum_dt < sim->m_max_t) )
        {            
            std::cout << "\n\n ------- sim step substrate ------- \n\n";
            std::cout << "curr t: " << sim->m_curr_t + accum_dt << " / " << sim->m_max_t << std::endl;
            
            std::cout<<"total vertices before improvement = "<<g_surf_subs->get_num_vertices()<<"\n";
            std::cout<<"total triangles before improvemen t= "<<g_surf_subs->m_mesh.num_triangles()<<"\n";
            
            // ---------- 
            // mesh maintenance & topology changes
            // ---------- 
            
            // Improve
            
            g_surf_subs->improve_mesh();
            g_surf_subs->defrag_mesh();
            
            std::cout<<"total vertices after improvement = "<<g_surf_subs->get_num_vertices()<<"\n";
            std::cout<<"total triangles after improvemen t= "<<g_surf_subs->m_mesh.num_triangles()<<"\n";            

            // ---------- 
            // advance underlying simulation
            // ----------
            
            // this is the desired time step for the sub-step, if that makes sense
            double curr_dt = sim_dt - accum_dt; 
            curr_dt = min( curr_dt, sim->m_max_t - sim->m_curr_t - accum_dt );
            
            std::cout << "curr_dt: " << curr_dt << std::endl;

            ////-------------------------------------------------------------////
            //*************** Explicit time integration ***********************//
            ////-------------------------------------------------------------////

            std::vector<Vec3d> new_positions( g_surf_subs->get_num_vertices() );
            // replace the driver with the explicit time integration with UL growth
            //driver->set_predicted_vertex_positions( *g_surf, new_positions, sim->m_curr_t + accum_dt, curr_dt );
            std::vector<Vec3d> V_t_hlf(g_surf_subs->get_num_vertices());
            //std::cout<<"explicit driver...\n";

            // July 13th 2019, off biaxial
            //std::cout<<"explicit driver...\n";
            //driverVelocities_Uniaxial(g_surf->get_positions(), g_surf->get_V_t(), V_t_hlf, g_surf->get_A_t(), curr_dt,sim->m_curr_t );
            //driverVelocities_OffBiaxial(g_surf->get_positions(), g_surf->get_V_t(), V_t_hlf, g_surf->get_A_t(), curr_dt,sim->m_curr_t );
            driverVelocities_Biaxial(g_surf_subs->get_positions(), g_surf_subs->get_V_t(), V_t_hlf, g_surf_subs->get_A_t(), curr_dt,sim->m_curr_t );
            //driverVelocities_BulgeTest(g_surf->get_positions(), g_surf->get_V_t(), V_t_hlf, g_surf->get_A_t(), curr_dt,sim->m_curr_t );
            //driverVelocities_Explicit(g_surf_subs->get_positions(), g_surf_subs->get_V_t(), V_t_hlf, g_surf_subs->get_A_t(), curr_dt );

            // use the V_t_hlf to get desired positions 
            for(int ni=0;ni<g_surf_subs->get_num_vertices();ni++){
                new_positions[ni] = g_surf_subs->get_position(ni)+curr_dt*V_t_hlf[ni];
                //std::cout<<"half velocity "<<V_t_hlf[ni]<<"\n";
                //std::cout<<" new["<<ni<<"]="<<new_positions[ni]<<", from old="<<g_surf->get_position(ni)<<"\n";                
            }
            //std::cout<<"done computing new positions based on V_hlf...\n";
            // this is setting the desired positions
            g_surf_subs->set_all_newpositions( new_positions );
            std::cout << "sim_dt: " << curr_dt << std::endl;
            
            // ----------
            // move & handle collision detection
            // ----------
            
            double actual_dt;
            std::vector<Vec3d> initial_positions = g_surf_subs->get_positions();
            
            // here the dynamic surface would like to indeed take the desired time step
            // but maybe it can't, and in that case it will return the actual time
            // step taken.
            std::cout<<"integrate...\n";
            g_surf_subs->integrate( curr_dt, actual_dt );
            
            // the time step actually taken by el topo
            double ori_dt = curr_dt;
            curr_dt = actual_dt;    
            
            std::cout << "actual_dt: " << actual_dt << std::endl;
            
            std::vector<Vec3d> final_positions = g_surf_subs->get_positions();
            //for(int ni=0;ni<final_positions.size();ni++){
            //    std::cout<<" final["<<ni<<"]="<<final_positions[ni]<<", from initial="<<initial_positions[ni]<<"\n";
            //}
            
            accum_dt += curr_dt;
            
            // with the actual final positions and initial positions recompute the residual
            // update the velocities, and also the strains
            // AND also pass the final vector of nodal reference integrin connections from the cell mesh
			
			// comment below for no_subs model
			
            updateForces_Explicit(initial_positions, final_positions, g_surf_subs->get_V_t(), V_t_hlf, g_surf_subs->get_A_t(), g_surf_subs->get_ea_t(), ori_dt, actual_dt, g_surf_subs->m_mesh, sim->m_curr_t, g_surf_cell->get_fa_u(),g_surf_cell->get_cC(),g_surf_cell->get_positions(),g_surf_cell->m_mesh);
            
            //
            // file output
            //
            
            //char binary_filename[256];
            //sprintf( binary_filename, "%s/frame%04d.bin", g_output_path, frame_stepper->get_frame() );      
            //write_binary_file( g_surf->m_mesh, g_surf->get_positions(), g_surf->m_masses, sim->m_curr_t, binary_filename );   

            //char vtk_filename[256];
            //sprintf( vtk_filename, "%s/subs_frame%04d.vtk", g_output_path, frame_stepper->get_frame() );      
            // several options
            //write_vtk( g_surf->m_mesh, g_surf->get_positions(), sim->m_curr_t, vtk_filename );   
            //write_vtk( g_surf_subs->m_mesh, g_surf_subs->get_positions(), g_surf_subs->get_cA(),g_surf_subs->get_cB(),g_surf_subs->get_cC(),sim->m_curr_t, vtk_filename );   
            //write_vtk( g_surf->m_mesh, g_surf->get_positions(), g_surf->get_ea_t(), sim->m_curr_t, vtk_filename );   
            
            //char stats_filename[256];
            //sprintf( stats_filename, "%s/aaa-stats.txt", g_output_path );      
            //g_stats.write_to_file( stats_filename );
            
        }
        
        std::cout << " -------------- end sim step substrate -------------- \n" << std::endl;

        // This function now is just the same as the matlab, takes in concentrations at 
        // previous time step and gives back new concentrations by reaction-diffusion 
        //advance_diffusion(g_surf->m_mesh, g_surf->get_positions(), g_surf->get_cB(),g_surf->get_cC(),g_surf->get_cBC(), sim_dt);

        
        sim->m_curr_t += accum_dt;
                
        
        // ----------
        // check if max time is reached
        // ----------
        
        if ( sim->done_simulation() )
        {         
            sim->m_running = false;
        //    std::cout << "total time steps: " << g_stats.get_int( "total_sim_steps" ) << std::endl;
        }
        
    }
    
    
    // ---------------------------------------------------------
    ///
    /// Advance simulation by one frame
    ///
    // ---------------------------------------------------------
    
    void advance_frame()
    {
        
        if( !sim->m_currently_advancing_simulation )
        {
            sim->m_currently_advancing_simulation = true;
            
            //
            // Advance frame
            //
            
            std::cout << "\n --------------- frame " << frame_stepper->get_frame() << " --------------- \n" << std::endl;
            
            while ( !(frame_stepper->done_frame()) )
            {
                double dt = frame_stepper->get_step_length( sim->m_dt );
                advance_sim( dt );
                // save at the end of the frame
                char vtk_cell_filename[256];
                sprintf( vtk_cell_filename, "%s/cell_frame%04d.vtk", g_output_path, frame_stepper->get_frame() );
				// double kint = -1; // use for MD-driven stiffness
				// double kint = 31426.8556; // use for linear stiffness of 32 pN/nm
				double kint = 1000; // use for 1 pN/nm
                write_vtk( g_surf_cell->m_mesh, g_surf_cell->get_positions(), g_surf_cell->get_cA(),g_surf_cell->get_cB(),g_surf_cell->get_cC(),g_surf_cell->get_fa_u(), g_surf_cell->get_ea_t(), sim->m_curr_t, vtk_cell_filename , kint);   
                char vtk_subs_filename[256];
                sprintf( vtk_subs_filename, "%s/subs_frame%04d.vtk", g_output_path, frame_stepper->get_frame() );      
                //write_vtk( g_surf_subs->m_mesh, g_surf_subs->get_positions(), g_surf_subs->get_cA(),g_surf_subs->get_cB(),g_surf_subs->get_cC(),sim->m_curr_t, vtk_subs_filename );   
                write_vtk( g_surf_subs->m_mesh, g_surf_subs->get_positions(), g_surf_subs->get_ea_t(), g_surf_subs->get_fa_u(), sim->m_curr_t, vtk_subs_filename );   
                frame_stepper->advance_step( dt ); 
                
            }
            
            std::cout << " --------------- end frame " << frame_stepper->get_frame() << " --------------- \n" << std::endl;
            
            // update frame stepper      
            frame_stepper->next_frame();
            
            sim->m_currently_advancing_simulation = false;

            // For some simulations, maybe there is need to update time step at different parts
            // of the simulation 
            //if(sim->m_curr_t>3){
            //    sim->m_dt = 0.0001;
            //}
            
        }
        
    }
    
    
    
    // ---------------------------------------------------------
    ///
    /// Run an entire simulation without GUI.  No threading.
    ///
    // ---------------------------------------------------------
    void run_simulation( )
    {

        sim->m_running = true;
        while ( sim->m_running )
        {
            advance_frame();
        }
        sim->m_running = false;
    }
    
    
    // ---------------------------------------------------------
    
    void set_up_output_path( )
    {
        
        DIR *dp;
        struct dirent *ep;     
        dp = opendir ( g_output_path );
        
        int max_dir = 0;
        
        if (dp != NULL)
        {      
            while ( (ep = readdir (dp)) != 0 )
            {
                
                std::string name_string( ep->d_name );
                
                //if ( ep->d_type == DT_DIR )    // doesn't exist on all platforms
                {
                    int i;
                    int num_variables_read = sscanf (name_string.c_str(), "run%04d", &i );
                    if ( num_variables_read > 0 )
                    {
                        max_dir = max( max_dir, i );
                    }
                }
            }
            
            (void) closedir (dp);
        }
        else
        {
            std::cout << "Couldn't open the specified output directory" << std::endl;
        }
        
        sprintf( g_output_path, "%s/run%04d", g_output_path, max_dir + 1 );
        std::cout << "Output path: " << g_output_path << std::endl;
        
        char mkdir_command[1024];
        sprintf( mkdir_command, "mkdir -p %s", g_output_path );
        system(mkdir_command);
                
    }
    
    
    // ---------------------------------------------------------
    ///
    /// Initialize the simulation.  Set the simulation parameters, initial geometry, etc.
    ///
    // ---------------------------------------------------------
    
    void init_simulation( char *script_filename_cell, char *script_filename_subs )
    {
        
        SurfTrackInitializationParameters init_params;
        
        // NEED to read two different surfaces. 
        // so two different script files 
        // First read the cell, with all the information relevant to the cell
        ScriptInit script_init;
        std::cout<<"reading cell script...\n";
        script_init.parse_script( script_filename_cell );
        std::cout<<"finished reading script\n";
        
        if ( script_init.output_path_is_relative )
        {
            snprintf( g_output_path, 256, "%s/%s", g_base_output_path, script_init.output_path.c_str() );
        }
        else
        {
            snprintf( g_output_path, 256, "%s", script_init.output_path.c_str() );
        }
        
        // Init frame stepper
        
        frame_stepper = new FrameStepper( script_init.frame_dt );
        
        // init simulation
        
        if ( script_init.end_sim_t != UNINITIALIZED_DOUBLE )
        {
            sim = new Simulation( script_init.sim_dt, script_init.end_sim_t );
        }
        else
        {
            sim = new Simulation( script_init.sim_dt );
        }
        
        if ( script_init.curr_t_specified )
        {
            sim->m_curr_t = script_init.curr_t;
            
            unsigned int curr_frame = static_cast<unsigned int>(script_init.curr_t / script_init.frame_dt); 
            std::cout << "curr_t: " << sim->m_curr_t << std::endl;
            std::cout << "curr_frame: " << curr_frame << std::endl; 
            
            frame_stepper->frame_count = curr_frame;
        }
        
        
        // init SurfTrack
        
        // NOTE: initializer for concentrations, velocities, accelerations 
        std::cout<<"initialize cell surfTrack\n";
        g_surf_cell = new SurfTrack( script_init.vertices, script_init.triangles, script_init.masses,script_init.c_A,script_init.c_B,script_init.c_C,script_init.V_t,script_init.A_t,script_init.ea_t,script_init.fa_u, script_init.surf_track_params );   
        std::cout<<"defrag mesh\n";
        g_surf_cell->defrag_mesh();

        // read another script but only for the mesh of the substrate, ignore the simulation parameters 
        ScriptInit script_init_subs;
        std::cout<<"reading subs script...\n";
        script_init_subs.parse_script( script_filename_subs );
        std::cout<<"finished reading script\n";
        
        
        std::cout<<"initialize subs surfTrack\n";
        g_surf_subs = new SurfTrack( script_init_subs.vertices, script_init_subs.triangles, script_init_subs.masses,script_init_subs.c_A,script_init_subs.c_B,script_init_subs.c_C,script_init_subs.V_t,script_init_subs.A_t,script_init_subs.ea_t,script_init.fa_u, script_init_subs.surf_track_params );   
        std::cout<<"defrag mesh\n";
        g_surf_subs->defrag_mesh();
        
        
        if ( g_surf_cell->m_collision_safety )
        {
            g_surf_cell->m_collision_pipeline.assert_mesh_is_intersection_free( false );      
        }

        if ( g_surf_subs->m_collision_safety )
        {
            g_surf_subs->m_collision_pipeline.assert_mesh_is_intersection_free( false );      
        }
        
        set_time_base();
        
    }
    
}  // unnamed namespace

// ---------------------------------------------------------
///
/// MAIN
///
// ---------------------------------------------------------

int main(int argc, char **argv)
{   
    
    std::cout << "\nCellBP: use of El Topo to simulate cell biophysics" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl << std::endl;
    
    if( argc < 2 )
    {
        std::cout << "Usage: <executable> <scriptfile_cell> <scriptfile_subs> <outputbasedirectory>" << std::endl;
        std::cout << "e.g.: " << std::endl;
        std::cout << "$ ./cellBP_release simulation_script_cell.txt simulation_script_subs.txt ./ \n" << std::endl;
        return 0;
    }
    
    // set path for outputting obj, bin files, etc.
    
    if ( argc > 3 )
    {
        strncpy( g_base_output_path, argv[3], sizeof(g_base_output_path) );
    }
    else
    {
        // argc == 3
        strncpy( g_base_output_path, "./", sizeof(g_base_output_path) );
    }
    

    //
    // Initialize the simulation using the script file
    // Basically read the mesh and the remeshing parameters for topo 
    //
    
    init_simulation( argv[1], argv[2] );
    
    
    //
    // Make a new directory for output
    //
    
    set_up_output_path();
    
    //
    // Make a copy of the input cell script in the output directory
    //
    
    char* script_filename_cell = argv[1];
    std::ifstream original_file_stream( script_filename_cell );
    assert( original_file_stream.good() );
    
    char script_copy_filename[256];
    snprintf( script_copy_filename, 256, "%s/aaa-script_cell.txt", g_output_path );
    
    char command[1024];
    sprintf( command, "cp %s %s", script_filename_cell, script_copy_filename );
    system(command);

    char* script_filename_subs = argv[2];
    std::ifstream original_file_stream2( script_filename_subs );
    assert( original_file_stream2.good() );
    
    char script_copy_filename2[256];
    snprintf( script_copy_filename2, 256, "%s/aaa-script_subs.txt", g_output_path );
    
    char command2[1024];
    sprintf( command2, "cp %s %s", script_filename_subs, script_copy_filename2 );
    system(command2);
    
    srand( 1 );
    
    // write frame 0
    //char binary_filename[256];
    //sprintf( binary_filename, "%s/frame%04d.bin", g_output_path, frame_stepper->get_frame() );      
    //write_binary_file( g_surf->m_mesh, g_surf->get_positions(), g_surf->m_masses, sim->m_curr_t, binary_filename );   
    
    char vtk_filename[256];
    sprintf( vtk_filename, "%s/cell_frame%04d.vtk", g_output_path, frame_stepper->get_frame() );      
    // double kint = -1; // use for MD-driven stiffness
	// double kint = 31426.8556; // use for linear stiffness
	double kint = 1000; // use for 1 pN/nm
	write_vtk( g_surf_cell->m_mesh, g_surf_cell->get_positions(), g_surf_cell->get_cA(),g_surf_cell->get_cB(),g_surf_cell->get_cC(), g_surf_cell->get_fa_u(), g_surf_cell->get_ea_t(), sim->m_curr_t, vtk_filename , kint);
    char vtk_filename2[256];
    sprintf( vtk_filename2, "%s/subs_frame%04d.vtk", g_output_path, frame_stepper->get_frame() );      
    //write_vtk( g_surf_subs->m_mesh, g_surf_subs->get_positions(), g_surf_subs->get_cA(),g_surf_subs->get_cB(),g_surf_subs->get_cC(), sim->m_curr_t, vtk_filename );
    write_vtk( g_surf_subs->m_mesh, g_surf_subs->get_positions(), g_surf_subs->get_ea_t(), g_surf_cell->get_fa_u(), sim->m_curr_t, vtk_filename2 );   
    //write_vtk( g_surf->m_mesh, g_surf->get_positions(), g_surf->get_ea_t(), sim->m_curr_t, vtk_filename );   

    //driver->write_to_disk( g_output_path, frame_stepper->get_frame() );
    
    //
    // Now start
    //
    
    // start the simulation (hands off)

    run_simulation();

    
    // uncomment to run all examples in our SISC paper:
    //run_all_sisc_examples( "./sisc-scripts/" );

    return 0;
}



