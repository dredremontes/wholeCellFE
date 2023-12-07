// ---------------------------------------------------------
//
//  scriptinit.h
//  Tyson Brochu 2011
//
//  Parse a script text file to initialize the simulation and mesh obects.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_SCRIPTINIT_H
#define EL_TOPO_SCRIPTINIT_H

#include <newparser.h>
#include <surftrack.h>
#include <vec.h>
#include <vector>

#include <Eigen/Dense>

class MeshDriver;

class ScriptInit
{
    
public:
    
    ScriptInit() :
    output_path(),
    output_path_is_relative(false),
    frame_dt( UNINITIALIZED_DOUBLE ),
    sim_dt( UNINITIALIZED_DOUBLE ),
    end_sim_t( UNINITIALIZED_DOUBLE ),
    curr_t_specified(false),
    curr_t(0.0),
    vertices(),
    triangles(),
    masses(),
    c_A(), c_B(), c_C(), V_t(),A_t(),ea_t(),fa_u(),a_i(),
    surf_track_params(),
    driver( NULL ),
    camera_target( UNINITIALIZED_DOUBLE ),
    camera_distance( UNINITIALIZED_DOUBLE ),
    camera_heading( UNINITIALIZED_DOUBLE ),
    camera_pitch( UNINITIALIZED_DOUBLE )
    {}
    
    void parse_script( const char* filename );
    
private:
    
    ScriptInit( const ScriptInit& other );
    const ScriptInit& operator=( const ScriptInit& other );
    
    void parse_surftrack_parameters( const ParseTree& surftrack_branch );
    
    // Mesh drivers
    // --------
    
    //void parse_invagination( const ParseTree& invagination_sim_branch );   

    void parse_camera( const ParseTree& camera_branch );
    
    // Geometry
    // --------
    
    void parse_sheet( const ParseTree& sheet_branch );
    void parse_curved_sheet( const ParseTree& curved_sheet_branch );
    void parse_sphere( const ParseTree& sphere_branch );
    void parse_dumbbell( const ParseTree& dumbbell_branch );
    
public:
    
    
    // Simulation settings
    // --------
    
    std::string output_path;
    bool output_path_is_relative;
    
    double frame_dt;
    double sim_dt;
    double end_sim_t;
    
    bool curr_t_specified;
    double curr_t;
        
    // Surface geometry
    // --------
    
    std::vector<Vec3d> vertices;
    std::vector<Vec3st> triangles;
    std::vector<double> masses;

    // NOTE: concentrations 
    std::vector<double> c_A;
    std::vector<double> c_B;
    std::vector<double> c_C;

    // NOTE: velocities, accelerations, strains
    std::vector<Vec3d> V_t;
    std::vector<Vec3d> A_t;
    std::vector<Eigen::VectorXd> ea_t;

    // NOTE: vector with the reference position of the substrate point connecting to FA
    std::vector<Vec3d> fa_u;

    std::vector<Vec3d> a_i; // Andre's a_i addition
    

    // SurfTrack
    // --------
    
    SurfTrackInitializationParameters surf_track_params;
    
    // MeshDriver
    // --------
    
    MeshDriver* driver;
    
    // GUI
    // --------
    
    Vec3d camera_target;
    double camera_distance, camera_heading, camera_pitch;
    
};

#endif


