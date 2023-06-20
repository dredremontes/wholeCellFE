#include "vec.h"
#include "nondestructivetrimesh.h"
#include <vector>
#include <math.h>   
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
//#include "densitytable.h"

//using namespace Eigen;

// Updated lagrangian formulation for growing membrane 

// function declarations
//void driverVelocities_Explicit(std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt);
//void updateForces_Explicit(std::vector<Vec3d> &X_t, std::vector<Vec3d> &x_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t, std::vector<Eigen::VectorXd> &ea_t, double dt, double dt_real, const NonDestructiveTriMesh &mesh);
void calculateRes(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Eigen::VectorXd> &ea_t,const std::vector<double> &c_A, const std::vector<double> &c_C, const NonDestructiveTriMesh &mesh, double dt_real, Eigen::VectorXd &Res, Eigen::VectorXd &Res_ea, Eigen::VectorXd &MM, Eigen::VectorXd &Res_Adif, Eigen::VectorXd &Res_Cdif, int mat, Eigen::Vector3d &centroid, double time);
void evalElementRe(Eigen::Vector3d &node1_X, Eigen::Vector3d &node2_X, Eigen::Vector3d &node3_X, Eigen::Vector3d &node1_x, Eigen::Vector3d &node2_x, Eigen::Vector3d &node3_x, Eigen::VectorXd &node1_ea, Eigen::VectorXd &node2_ea, Eigen::VectorXd &node3_ea, std::vector<double> &node_c_A, std::vector<double> &node_c_C, double dt_real, Eigen::VectorXd &Re, Eigen::VectorXd &Re_ea, Eigen::VectorXd &Re_Adif, Eigen::VectorXd &Re_Cdif, Eigen::Vector3d &MMe, int mat, Eigen::Vector3d &centroid, double time);
void calculateRes(const std::vector<Vec3d> &X_t,const std::vector<Vec3d> &x_t, std::vector<Eigen::VectorXd> &ea_t,const NonDestructiveTriMesh &mesh, double dt_real, Eigen::VectorXd &Res, Eigen::VectorXd &Res_ea, Eigen::VectorXd &MM, int mat, double time);
void evalElementRe(Eigen::Vector3d &node1_X, Eigen::Vector3d &node2_X, Eigen::Vector3d &node3_X, Eigen::Vector3d &node1_x, Eigen::Vector3d &node2_x, Eigen::Vector3d &node3_x, Eigen::VectorXd &node1_ea, Eigen::VectorXd &node2_ea, Eigen::VectorXd &node3_ea, double dt_real, Eigen::VectorXd &Re, Eigen::VectorXd &Re_ea, Eigen::Vector3d &MMe, int mat, double time);
void calculateConstraintForce(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &x_t);
void calculatePressureForce(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &x_t, const Eigen::VectorXd &MM, const NonDestructiveTriMesh &mesh, double t_tot);
void calculateBendingRes(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const NonDestructiveTriMesh &mesh, Eigen::VectorXd &Res_bending);
double calculateBendingEnergy(const std::vector<Eigen::Vector3d> &v_vec_X , const std::vector<Eigen::Vector3d> &v_vec_x,bool print_info);
void calculateForce_FA(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const std::vector<double> &c_C, std::vector<Vec3d> &fa_u, const NonDestructiveTriMesh &mesh, double kint);
void calculateForce_FA_subs(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const NonDestructiveTriMesh &mesh, const std::vector<double> &c_C, const std::vector<Vec3d> &x_cell, std::vector<Vec3d> &fa_u, const NonDestructiveTriMesh &mesh_cell,double kint);
void calculateForce_AP(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const std::vector<double> &c_A, const NonDestructiveTriMesh &mesh, double dt_real);
void calculateBendingForce2D(Eigen::VectorXd &F_bending, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const NonDestructiveTriMesh &mesh);
Eigen::Vector4d calculateAreaTot(const std::vector<Vec3d> &x_t,const NonDestructiveTriMesh &mesh);
void calculateAreaConstraintForce(Eigen::VectorXd &F_areaConstraint, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const NonDestructiveTriMesh &mesh, double tot_area);
void updateFA(std::vector<double> &c_C, std::vector<Vec3d> &fa_u, double dt_real, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t,double kint);
void updateFA(std::vector<double> &c_A, std::vector<double> &c_C, std::vector<Vec3d> &fa_u, double dt_real, const std::vector<Vec3d> &X_t_cell, const std::vector<Vec3d> &x_t_cell,const NonDestructiveTriMesh &mesh_cell, const std::vector<Vec3d> &X_t_subs, const std::vector<Vec3d> &x_t_subs,const NonDestructiveTriMesh &mesh_subs, double kint);

// Driver for Zebrafish based on reading from file
void driverVelocities_Zebrafish(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt, double t_tot){
	// Open the corresponding file with the velocities and read into vectors
	int currframe = 4*int(t_tot/6)+130; // each frame is 1.5 min 
	std::cout<<"frame "<<currframe<<"\n";
	std::string filename="velmap_reflexed/velocity_at";
	filename.append(std::to_string(currframe));
	filename.append(".txt");
	std::ifstream velocityfile(filename.c_str());
	std::cout<<"file: "<<filename<<"\n";
    std::string line;
    std::vector<Vec3d> cell_positions; cell_positions.clear(); // positions in microns 
    std::vector<Vec3d> cell_u; cell_u.clear(); // displacement in microns 
    if (velocityfile.is_open())
    {
        while ( getline (velocityfile,line) )
        {
        	std::stringstream ss0(line);
        	std::vector<std::string> result;
			while( ss0.good() )
			{
    			std::string substr;
    			getline( ss0, substr, ',' );
    			result.push_back( substr );
			}
            Vec3d new_cell, new_u; 
            new_cell[0]= std::stod(result[0]);
            new_cell[1]= std::stod(result[1]);
            new_cell[2]= std::stod(result[2]);
            new_u[0]= std::stod(result[3]);
            new_u[1]= std::stod(result[4]);
            new_u[2]= std::stod(result[5]);
            cell_positions.push_back(new_cell);
            cell_u.push_back(new_u);
            //std::cout<<new_cell<<"\n";
            //std::cout<<new_u<<"\n\n";
        }
    }
    // Loop over the nodes of the mesh and find the closest cells and compute their average velocity
    int ncell = cell_positions.size();
    double radius = 15; // microns
    std::vector<int> neighbors;neighbors.clear();
    int nn = X_t.size();
	for(int ni=0;ni<nn;ni++){
		for(int nj=0;nj<ncell;nj++){
			Vec3d d_vec = X_t[ni]-cell_positions[nj];
			double d_mag = sqrt(d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]);
			if(d_mag<radius){
				neighbors.push_back(nj);
			}
		}
		// given the neighbor list get their average velocity 
		Vec3d vel_aux; vel_aux[0]=0.;vel_aux[1]=0.;vel_aux[2]=0.;
		for(int nj=0;nj<neighbors.size();nj++){
			vel_aux = vel_aux + cell_u[neighbors[nj]]/6.;
			// note: displacements are in microns every 6min 
			// velocity is in um/min
		}
		if(neighbors.size()>1){vel_aux = vel_aux/neighbors.size();}

		// Update the velocities of the node 
		V_t_hlf[ni] = vel_aux; // um/min
		//std::cout<<vel_aux<<"\n";
	}
}

void driverVelocities_Explicit(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt){

	// Central difference following Belytshko Page 333

	// X_t contains the positions at time t
	// V_t contains the velocities at time t
	// A_t contains the velocities at time t
	// ea_t contains the euler almansi strain per node at time t (not needed in this function)
	// dt is the desired time step, note that Topo might take smaller time step

	// Velocities that we are going to give to Topo
	// These are velocities at the half time step
	int nn = X_t.size();
	std::cout<<"loop over nodes to get half velocities...\n";
	for(int ni=0;ni<nn;ni++){
		V_t_hlf[ni] = V_t[ni] + (dt/2.0)*A_t[ni] ;
		// turn off the normal component for the planar examples
		V_t_hlf[ni][2]=0.0;
	}
	
	std::cout<<"done.\n";
}
void driverVelocities_Uniaxial(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt, double t_tot){

	// Central difference following Belytshko Page 333

	// X_t contains the positions at time t
	// V_t contains the velocities at time t
	// A_t contains the velocities at time t
	// ea_t contains the euler almansi strain per node at time t (not needed in this function)
	// dt is the desired time step, note that Topo might take smaller time step

	// Velocities that we are going to give to Topo
	// These are velocities at the half time step
	int nn = X_t.size();
	std::cout<<"loop over nodes to get half velocities...\n";
	int n_cont_v = 0;
	int n_cont_s = 0;
	for(int ni=0;ni<nn;ni++){
		V_t_hlf[ni] = V_t[ni] + (dt/2.0)*A_t[ni] ;
		// turn off the normal component for the planar example
		V_t_hlf[ni][2]=0.0;
		// For the left fix and for the right impose a deformation 
		if(X_t[ni][0]<-0.0499){
			//std::cout<<"constrain left.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][0]>0.9499+0.1*t_tot && t_tot<=2.0){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.1;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
			n_cont_v+=1;
		}
		if(X_t[ni][0]>0.9499+2.00*0.1 && t_tot>=2.0){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
			n_cont_s+=1;
		}
	}
	
	std::cout<<"done, imposed "<<n_cont_v<<", "<<n_cont_s<<" right constraints.\n";
}
///
void driverVelocities_BulgeTest(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt, double t_tot){

	// Central difference following Belytshko Page 333

	// X_t contains the positions at time t
	// V_t contains the velocities at time t
	// A_t contains the velocities at time t
	// ea_t contains the euler almansi strain per node at time t (not needed in this function)
	// dt is the desired time step, note that Topo might take smaller time step

	// Velocities that we are going to give to Topo
	// These are velocities at the half time step
	int nn = X_t.size();
	std::cout<<"loop over nodes to get half velocities...\n";

	for(int ni=0;ni<nn;ni++){
		V_t_hlf[ni] = V_t[ni] + (dt/2.0)*A_t[ni] ;
		// turn off the normal component for the planar example
		//V_t_hlf[ni][2]=0.0;
		// For the left fix and for the right impose a deformation 
		if(X_t[ni][0]<0.001 && X_t[ni][2]<0.001 && X_t[ni][2]>-0.001){
			//std::cout<<"constrain left.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][0]>0.999 && X_t[ni][2]<0.001 && X_t[ni][2]>-0.001){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][1]<0.001 && X_t[ni][2]<0.001 && X_t[ni][2]>-0.001){
			//std::cout<<"constrain bottom.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][1]>0.999 && X_t[ni][2]<0.001 && X_t[ni][2]>-0.001){
			//std::cout<<"constrain top.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
	}
	
	//std::cout<<"done, imposed "<<n_cont_v<<", "<<n_cont_s<<" right constraints.\n";
}
/// 
void driverVelocities_OffBiaxial(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt, double t_tot){

	// Central difference following Belytshko Page 333

	// X_t contains the positions at time t
	// V_t contains the velocities at time t
	// A_t contains the velocities at time t
	// ea_t contains the euler almansi strain per node at time t (not needed in this function)
	// dt is the desired time step, note that Topo might take smaller time step

	// Velocities that we are going to give to Topo
	// These are velocities at the half time step
	int nn = X_t.size();
	std::cout<<"loop over nodes to get half velocities...\n";
	int n_cont_v = 0;
	int n_cont_s = 0;
	for(int ni=0;ni<nn;ni++){
		V_t_hlf[ni] = V_t[ni] + (dt/2.0)*A_t[ni] ;
		// turn off the normal component for the planar example
		V_t_hlf[ni][2]=0.0;
		// For the left fix and for the right impose a deformation 
		if(X_t[ni][0]<-0.04){
			//std::cout<<"constrain left.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][0]>0.94+0.1*t_tot && t_tot<1.0){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.1;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
			n_cont_v+=1;
		}
		if(X_t[ni][0]>0.94+0.99*0.1 && t_tot>=1.0){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
			n_cont_s+=1;
		}
		if(X_t[ni][1]>0.99){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][1]<0.01){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
	}
	
	//std::cout<<"done, imposed "<<n_cont_v<<", "<<n_cont_s<<" right constraints.\n";
}
void driverVelocities_Biaxial(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt, double t_tot){

	// Central difference following Belytshko Page 333

	// X_t contains the positions at time t
	// V_t contains the velocities at time t
	// A_t contains the velocities at time t
	// ea_t contains the euler almansi strain per node at time t (not needed in this function)
	// dt is the desired time step, note that Topo might take smaller time step

	// Velocities that we are going to give to Topo
	// These are velocities at the half time step
	int nn = X_t.size();
	std::cout<<"loop over nodes to get half velocities...\n";
	int n_cont_v = 0;
	int n_cont_s = 0;
	double Vbiax_x = 0.0; //um/s, initial geometry is 100um, so 1um in 1 s is 1%/sec, or 0.01 sec^-1 strain, ok
	double Vbiax_y = 0.0; 
	// only do the biaxial after some equilibration of 3 sec
	for(int ni=0;ni<nn;ni++){
		V_t_hlf[ni] = V_t[ni] + (dt/2.0)*A_t[ni] ;
		// turn off the normal component for the planar example
		V_t_hlf[ni][2]=0.0;
		double t_i = 3;
		double t_f = 8;
		if(t_tot<t_i){
			if(X_t[ni][0]<-49.00){
				V_t_hlf[ni][0]=0.0;
			}
			if(X_t[ni][0]>49.00){
				V_t_hlf[ni][0]=0.0;
			}
			if(X_t[ni][1]>49.00){
				V_t_hlf[ni][1]=0.0;
			}
			if(X_t[ni][1]<-49.00){
				V_t_hlf[ni][1]=0.0;
			}
		}
		if(t_tot>t_i && t_tot<t_f){
			if(X_t[ni][0]<-49.00-(t_tot-t_i)*Vbiax_x){
				V_t_hlf[ni][0]=-Vbiax_x;
			}
			if(X_t[ni][0]>49.00+(t_tot-t_i)*Vbiax_x){
				V_t_hlf[ni][0]=Vbiax_x;
			}
			if(X_t[ni][1]>49.00+(t_tot-t_i)*Vbiax_y){
				V_t_hlf[ni][1]=Vbiax_y;
			}
			if(X_t[ni][1]<-49.00-(t_tot-t_i)*Vbiax_y){
				V_t_hlf[ni][1]=-Vbiax_y;
			}
		}
		if(t_tot>t_f){
			if(X_t[ni][0]<-49.00-(t_f-t_i)*Vbiax_x){
				V_t_hlf[ni][0]=0.0;
			}
			if(X_t[ni][0]>49.00+(t_f-t_i)*Vbiax_x){
				V_t_hlf[ni][0]=0.0;
			}
			if(X_t[ni][1]>49.00+(t_f-t_i)*Vbiax_y){
				V_t_hlf[ni][1]=0.0;
			}
			if(X_t[ni][1]<-49.00-(t_f-t_i)*Vbiax_y){
				V_t_hlf[ni][1]=0.0;
			}
		}
	}
	
	//std::cout<<"done, imposed "<<n_cont_v<<", "<<n_cont_s<<" right constraints.\n";
}

// This is the default function, doesnt include cell-substrate interactions or anything like that
void updateForces_Explicit(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t, std::vector<Eigen::VectorXd> &ea_t, double dt, double dt_real, const NonDestructiveTriMesh &mesh, double t_tot){

	// Central difference following Belytshko Page 333 

	// at this point we passed the velocities of the mid time point to Topo
	// and Topo integrated the equation and gives us the actual current positions
	// and the actual time step, same mesh as before, the mesh changes are done
	// separately

	//** Update the force f_t  **//
	// calculate the internal residual
	int n_node = X_t.size();
	Eigen::VectorXd Res(n_node*3);
	// I also need the residual for the update of the strain
	Eigen::VectorXd Res_ea(n_node*6);
	// And here I can get the lumped mass 'matrix' easily as well
	Eigen::VectorXd MM(n_node);

	std::cout<<"going into calculateRes() \n";
	std::cout<<"number of nodes= "<<n_node<<"\n";
	//******************************************************//
	// passing 'subs', but this mostly passed because needs to be passed
	// however, as mentioned at the start of the function, this function was the original one
	// doesnt do anything in terms of the area constraint, bending energy, FA forces, etc.
	// but still need to pass which 'material', so default substrate
	int mat = 2; 
	calculateRes(X_t,x_t,ea_t,mesh,dt_real, Res,Res_ea,MM,mat,t_tot); 
	//******************************************************//

	// Calculate the force constraint, forces needed to enforce things like:
	// * growing tissue wants to stick to the sphere, or similar forces
	// Eigen::VectorXd F_constraint(n_node*3);F_constraint.setZero();
	//calculateConstraintForce(F_constraint,x_t);
	
	// Pressure force for the bulge test example
	//calculatePressureForce(F_constraint,x_t,MM, mesh, t_tot);

	// Calculate bending residual, to prevent buckling, treat as regularization
	//Eigen::VectorXd Res_bending(n_node*3);Res_bending.setZero();
	//calculateBendingRes(X_t,x_t,mesh,Res_bending);
	//std::cout<<"bending residual\n"<<Res_bending<<"\n";

	// net force
	//Eigen::VectorXd f_t = F_constraint - Res - Res_bending;
	Eigen::VectorXd f_t = - Res;
	std::cout<<"net force norm "<<Res.norm()<<"\n";
	// Update the Acceleration for the next time step, and also the strain
	double damp = 2.0; // [(kg/s)] not sure this is good or bad? 1 was ok for uniaxial, 0.1 was ok for ZF
	//double rho = 1.0; // [ug/um^3] also not sure if this is good or not
	double rho = 1000.0; // [kg/m^3] 
	// Update acceleration and strain
	for(int ni=0;ni<n_node;ni++){
		// updating the acceleration
		A_t[ni][0] = (f_t(ni*3+0)-damp*MM(ni)*rho*V_t[ni][0])/(MM(ni)*rho);
		A_t[ni][1] = (f_t(ni*3+1)-damp*MM(ni)*rho*V_t[ni][1])/(MM(ni)*rho);
		A_t[ni][2] = (f_t(ni*3+2)-damp*MM(ni)*rho*V_t[ni][2])/(MM(ni)*rho);	
		// updating the strains at the nodes
		ea_t[ni](0) = Res_ea(ni*6+0)/MM(ni);
		ea_t[ni](1) = Res_ea(ni*6+1)/MM(ni);
		ea_t[ni](2) = Res_ea(ni*6+2)/MM(ni);
		ea_t[ni](3) = Res_ea(ni*6+3)/MM(ni);
		ea_t[ni](4) = Res_ea(ni*6+4)/MM(ni);
		ea_t[ni](5) = Res_ea(ni*6+5)/MM(ni);
		//std::cout<<"node "<<ni<<", f_t = "<<f_t(ni*3+0)<<", "<<f_t(ni*3+1)<<", "<<f_t(ni*3+2)<<"\n";
		//std::cout<<"node "<<ni<<", Res = "<<Res(ni*3+0)<<", "<<Res(ni*3+1)<<", "<<Res(ni*3+2)<<"\n";
	}
	
	// Update the velocities for next time step
	for(int ni=0;ni<n_node;ni++){
		if(dt_real>dt/2.0){
			// this should always be the case in that most of the time dt_real = dt,
			// we can print is dt_real and dt are not the same... 
			V_t[ni] = V_t_hlf[ni] + (dt_real - dt/2.0)*A_t[ni];
		}else{
			// if dt_real < dt, then this might stil be too large 
			V_t[ni] = V_t[ni] + dt_real*A_t[ni];
		}
	}
}

// This function includes the vector of focal adhesion attachements to model cell-subs interaction
void updateForces_Explicit(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t, std::vector<Eigen::VectorXd> &ea_t, double dt, double dt_real, const NonDestructiveTriMesh &mesh, double t_tot, std::vector<Vec3d> &fa_u, const std::vector<double> &c_C, const std::vector<Vec3d> &x_cell, const NonDestructiveTriMesh &mesh_cell){

	// Central difference following Belytshko Page 333 

	// at this point we passed the velocities of the mid time point to Topo
	// and Topo integrated the equation and gives us the actual current positions
	// and the actual time step, same mesh as before, the mesh changes are done
	// separately

	//** Update the force f_t  **//
	// calculate the internal residual
	int n_node = X_t.size();
	Eigen::VectorXd Res(n_node*3);
	// I also need the residual for the update of the strain
	Eigen::VectorXd Res_ea(n_node*6);
	// And here I can get the lumped mass 'matrix' easily as well
	Eigen::VectorXd MM(n_node);

	std::cout<<"going into calculateRes() for substrate \n";
	std::cout<<"number of nodes= "<<n_node<<"\n";

	//******************************************************//
	int mat = 2; // for substrate 
	calculateRes(X_t,x_t,ea_t,mesh,dt_real, Res,Res_ea,MM,mat,t_tot); 
	//******************************************************//

	// Calculate the force from the cells on the substrate
	Eigen::VectorXd F_FA(n_node*3);F_FA.setZero();
	double kint;
	kint = 1000; // 1000 pN/um = 1 pN/nm
	// kint = 31000;
	// if(t_tot<4){
	//	kint = 1000; // 1000 pN/um = 1 pN/nm 
	//}else if(t_tot<8){
	//	kint = 10000;
	//}else if(t_tot<12){
	//	kint = 100000;
	// }
	calculateForce_FA_subs(F_FA,X_t,x_t,mesh,c_C,x_cell,fa_u,mesh_cell,kint);

	// net force
	//Eigen::VectorXd f_t = F_constraint - Res - Res_bending;
	Eigen::VectorXd f_t = F_FA - Res;
	std::cout<<"net internal force norm on substrate "<<Res.norm()<<"\n";
	// Update the Acceleration for the next time step, and also the strain
	double damp = 0.001; // [(1/s)] not sure this is good or bad? 1 was ok for uniaxial, 0.1 was ok for ZF
	double rho = 1.0; // [ug/um^3] also not sure if this is good or not
	//double rho = 1.0; // [pg/um^3] 
	// Update acceleration and strain
	for(int ni=0;ni<n_node;ni++){
		// updating the acceleration
		A_t[ni][0] = (f_t(ni*3+0)-damp*MM(ni)*rho*V_t[ni][0])/(MM(ni)*rho);
		A_t[ni][1] = (f_t(ni*3+1)-damp*MM(ni)*rho*V_t[ni][1])/(MM(ni)*rho);
		A_t[ni][2] = (f_t(ni*3+2)-damp*MM(ni)*rho*V_t[ni][2])/(MM(ni)*rho);	
		// updating the strains at the nodes
		// June 2, 2023, change to a different update based on incremental strain rather than 
		// complete removal 
		ea_t[ni](0) += Res_ea(ni*6+0)/MM(ni);
		ea_t[ni](1) += Res_ea(ni*6+1)/MM(ni);
		ea_t[ni](2) += Res_ea(ni*6+2)/MM(ni);
		ea_t[ni](3) += Res_ea(ni*6+3)/MM(ni);
		ea_t[ni](4) += Res_ea(ni*6+4)/MM(ni);
		ea_t[ni](5) += Res_ea(ni*6+5)/MM(ni);
		//std::cout<<"node "<<ni<<", f_t = "<<f_t(ni*3+0)<<", "<<f_t(ni*3+1)<<", "<<f_t(ni*3+2)<<"\n";
		//std::cout<<"node "<<ni<<", Res = "<<Res(ni*3+0)<<", "<<Res(ni*3+1)<<", "<<Res(ni*3+2)<<"\n";
	}
	
	// Update the velocities for next time step
	for(int ni=0;ni<n_node;ni++){
		if(dt_real>dt/2.0){
			V_t[ni] = V_t_hlf[ni] + (dt_real - dt/2.0)*A_t[ni];
		}else{
			V_t[ni] = V_t[ni] + dt_real*A_t[ni];
		}
	}
}


// This is the overwritten version for the cell because it includes volume control, bending of membrane
void updateForces_Explicit(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t, std::vector<Eigen::VectorXd> &ea_t, 
	std::vector<double> &c_A, std::vector<double> &c_B, std::vector<double> &c_C, std::vector<Vec3d> &fa_u, double dt, double dt_real, const NonDestructiveTriMesh &mesh, double t_tot){

	// Central difference following Belytshko Page 333 

	// at this point we passed the velocities of the mid time point to Topo
	// and Topo integrated the equation and gives us the actual current positions
	// and the actual time step, same mesh as before, the mesh changes are done
	// separately

	//** Update the force f_t  **//
	// calculate the internal residual
	int n_node = X_t.size();
	Eigen::VectorXd Res(n_node*3);
	// I also need the residual for the update of the strain
	Eigen::VectorXd Res_ea(n_node*6);
	// And here I can get the lumped mass 'matrix' easily as well
	Eigen::VectorXd MM(n_node);
	// For diffusion 
	Eigen::VectorXd Res_Adif(n_node);
	Eigen::VectorXd Res_Cdif(n_node);

	// before going into residual function, for anisotropic stress 
	Eigen::Vector4d area = calculateAreaTot(x_t,mesh);
	double area_tot = area(0);
	Eigen::Vector3d centroid; centroid<<area(1),area(2),area(3);
	std::cout<<"total area "<<area_tot<<"\n";
	std::cout<<"centroid "<<centroid(0)<<", "<<centroid(1)<<", "<<centroid(2)<<"\n";

	std::cout<<"going into calculateRes() for cell \n";
	std::cout<<"number of nodes= "<<n_node<<"\n";

	// This residual includes both passive and active stresses 
	//******************************************************//
	int mat = 1; // for cell
	calculateRes(X_t,x_t,ea_t, c_A, c_C, mesh,dt_real, Res,Res_ea, MM, Res_Adif, Res_Cdif, mat, centroid, t_tot); 
	//******************************************************//

	// Calculate the force constraint, forces needed to enforce things like:
	// * growing tissue wants to stick to the sphere
	Eigen::VectorXd F_FA(n_node*3);F_FA.setZero();
	Eigen::VectorXd F_AP(n_node*3);F_AP.setZero();
	
	//calculateConstraintForce(F_constraint,x_t);

	// Force for focal adhesions, need previous and new positions
	// as well as whether there is a FA, maybe will need additional parameters
	// or state variables for each adhesion, but first just one
	// Update focal adhesion concentration
	std::cout<<"going to eval constraint from FA\n";
	double kint;
	kint = 1000; // 1000 pN/um = 1 pN/nm
	// kint = 31000;

	calculateForce_FA(F_FA,X_t,x_t,c_C,fa_u,mesh,kint);

	// Also need to add some actin polymerization force 
	calculateForce_AP(F_AP,X_t,x_t,c_A,mesh,dt_real);	

	// add some viscous drag for stability 
	Eigen::VectorXd F_vd(n_node*3);F_vd.setZero();
	double drag = 0.001; // pN*s/um 
	for(int ni=0;ni<x_t.size();ni++){
		F_vd(ni*3+0) = -drag*(V_t_hlf[ni][0]);
		F_vd(ni*3+1) = -drag*(V_t_hlf[ni][1]);
		F_vd(ni*3+2) = -drag*(V_t_hlf[ni][2]);
	}

	// Pressure force for the bulge test
	//calculatePressureForce(F_constraint,x_t,MM, mesh, t_tot);

	// Calculate bending residual, to prevent buckling, treat as regularization
	Eigen::VectorXd F_bending(n_node*3);F_bending.setZero();
	//calculateBendingRes(X_t,x_t,mesh,Res_bending);
	// this specifically is for the curvature of a curve embedded in 2D, i.e. the curvature of
	// the boundary of the mesh 
	calculateBendingForce2D(F_bending,X_t,x_t,mesh);
	//std::cout<<"bending residual\n"<<Res_bending<<"\n";

	// Area constraint, pressure-like force with respect to total area
	Eigen::VectorXd F_areaConstraint(n_node*3); F_areaConstraint.setZero(); 
	calculateAreaConstraintForce(F_areaConstraint, X_t, x_t, mesh, area_tot);

	// net force
	Eigen::VectorXd f_t = F_vd + F_bending + F_AP + F_FA + F_areaConstraint - Res;
	std::cout<<"net internal force norm "<<Res.norm()<<"\n";
	// Update the Acceleration for the next time step, and also the strain
	//double damp = 2.0; // [(kg/s)] not sure this is good or bad? 1 was ok for uniaxial, 0.1 was ok for ZF
	double damp = 0.001; // 1/s
	double rho = 1.0; // [ug/um^3] also not sure if this is good or not
	//double rho = 1000.0; // [kg/m^3] 
	//double rho = 1.0; // [pg/um^3] 
	// Update acceleration and strain
	for(int ni=0;ni<n_node;ni++){
		// updating the acceleration
		A_t[ni][0] = (f_t(ni*3+0)-damp*MM(ni)*rho*V_t[ni][0])/(MM(ni)*rho);
		A_t[ni][1] = (f_t(ni*3+1)-damp*MM(ni)*rho*V_t[ni][1])/(MM(ni)*rho);
		A_t[ni][2] = (f_t(ni*3+2)-damp*MM(ni)*rho*V_t[ni][2])/(MM(ni)*rho);	
		// updating the strains at the nodes
		ea_t[ni](0) += Res_ea(ni*6+0)/MM(ni);
		ea_t[ni](1) += Res_ea(ni*6+1)/MM(ni);
		ea_t[ni](2) += Res_ea(ni*6+2)/MM(ni);
		ea_t[ni](3) += Res_ea(ni*6+3)/MM(ni);
		ea_t[ni](4) += Res_ea(ni*6+4)/MM(ni);
		ea_t[ni](5) += Res_ea(ni*6+5)/MM(ni);
		//std::cout<<"node "<<ni<<", f_t = "<<f_t(ni*3+0)<<", "<<f_t(ni*3+1)<<", "<<f_t(ni*3+2)<<"\n";
		//std::cout<<"node "<<ni<<", Res = "<<Res(ni*3+0)<<", "<<Res(ni*3+1)<<", "<<Res(ni*3+2)<<"\n";
		// 
		// c_A, c_C
		//c_A[ni] = Res_Adif(ni)/MM(ni);
		//c_C[ni] = Res_Cdif(ni)/MM(ni);
	}
	
	// Update the velocities for next time step
	for(int ni=0;ni<n_node;ni++){
		if(dt_real>dt/2.0){
			V_t[ni] = V_t_hlf[ni] + (dt_real - dt/2.0)*A_t[ni];
		}else{
			V_t[ni] = V_t[ni] + dt_real*A_t[ni];
		}
	}

}

void calculateRes(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Eigen::VectorXd> &ea_t,const std::vector<double> &c_A, const std::vector<double> &c_C, const NonDestructiveTriMesh &mesh, double dt_real, Eigen::VectorXd &Res, Eigen::VectorXd &Res_ea, Eigen::VectorXd &MM, Eigen::VectorXd &Res_Adif, Eigen::VectorXd &Res_Cdif, int mat, Eigen::Vector3d &centroid, double time){
	Res.setZero();
	Res_ea.setZero();
	MM.setZero();
	// June 2, 2023
	// ABT adding diffusion, first easiest to do is to just diffuse c_C with small diffusion  
	// but then we actually need to also add reactions to transform between c_A and c_C
	Res_Adif.setZero();
	Res_Cdif.setZero();
	// START LOOP OVER ELEMENTS
	Eigen::VectorXd Re(9);
	Eigen::VectorXd Re_Adif(3);
	Eigen::VectorXd Re_Cdif(3);
	Eigen::VectorXd Re_ea(18);
	Eigen::Vector3d MMe;
    int n_elem = mesh.num_triangles();
    for(int ei=0;ei<n_elem;ei++){
       	// connectivity of the linear elements
        int node1index = mesh.get_triangle(ei)[0];
    	int node2index = mesh.get_triangle(ei)[1];
    	int node3index = mesh.get_triangle(ei)[2];
        
        // previous time step nodal positions for this element
		Eigen::Vector3d node1_X; node1_X<<X_t[node1index][0],X_t[node1index][1],X_t[node1index][2];
		Eigen::Vector3d node2_X; node2_X<<X_t[node2index][0],X_t[node2index][1],X_t[node2index][2];
		Eigen::Vector3d node3_X; node3_X<<X_t[node3index][0],X_t[node3index][1],X_t[node3index][2];

		// previous time step euler almansi strains at the nodes
		Eigen::VectorXd node1_ea(6); 
		node1_ea<<ea_t[node1index](0),ea_t[node1index](1),ea_t[node1index](2),ea_t[node1index](3),ea_t[node1index](4),ea_t[node1index](5);
		Eigen::VectorXd node2_ea(6); 
		node2_ea<<ea_t[node2index](0),ea_t[node2index](1),ea_t[node2index](2),ea_t[node2index](3),ea_t[node2index](4),ea_t[node2index](5);
		Eigen::VectorXd node3_ea(6); 
		node3_ea<<ea_t[node3index](0),ea_t[node3index](1),ea_t[node3index](2),ea_t[node3index](3),ea_t[node3index](4),ea_t[node3index](5);

		// current nodal positions for this element
		Eigen::Vector3d node1_x; node1_x<<x_t[node1index][0],x_t[node1index][1],x_t[node1index][2];
		Eigen::Vector3d node2_x; node2_x<<x_t[node2index][0],x_t[node2index][1],x_t[node2index][2];
		Eigen::Vector3d node3_x; node3_x<<x_t[node3index][0],x_t[node3index][1],x_t[node3index][2];

		// concentrations 
		std::vector<double> node_c_A(3); 
		node_c_A[0]=c_A[node1index]; node_c_A[1]=c_A[node2index]; node_c_A[2]=c_A[node3index];
		std::vector<double> node_c_C(3); 
		node_c_C[0]=c_C[node1index]; node_c_C[1]=c_C[node2index]; node_c_C[2]=c_C[node3index];
						
        // and calculate the element Re
        Re.setZero();
        Re_ea.setZero();
        MMe.setZero();
        // element residual for diffusion 
        Re_Adif.setZero();
        Re_Cdif.setZero();

				
        // subroutine to evaluate the element
        //std::cout<<"evaluating element "<<ei<<"\n";
        evalElementRe(node1_X,node2_X,node3_X,node1_x,node2_x,node3_x,node1_ea,node2_ea,node3_ea, node_c_A, node_c_C, dt_real, Re, Re_ea, Re_Adif, Re_Cdif, MMe, mat, centroid, time);
        //std::cout<<"element residual norm "<<Re.norm()<<"\n";
        // LOOP OVER NODES
		for(int nodei=0;nodei<3;nodei++){
			// ASSEMBLE DISPLACEMENT RESIDUAL
			for(int coordi=0;coordi<3;coordi++){
				// residual
				Res(mesh.get_triangle(ei)[nodei]*3+coordi) += Re(nodei*3+coordi);
			}
			for(int coordi=0;coordi<6;coordi++){
				Res_ea(mesh.get_triangle(ei)[nodei]*6+coordi) += Re_ea(nodei*6+coordi);
			}
			MM(mesh.get_triangle(ei)[nodei]) += MMe(nodei);
			// diffusion 
			Res_Adif(mesh.get_triangle(ei)[nodei]) += Re_Adif(nodei);
			Res_Cdif(mesh.get_triangle(ei)[nodei]) += Re_Cdif(nodei);
		}
	}
}


void evalElementRe(Eigen::Vector3d &node1_X, Eigen::Vector3d &node2_X, Eigen::Vector3d &node3_X, Eigen::Vector3d &node1_x, Eigen::Vector3d &node2_x, Eigen::Vector3d &node3_x, Eigen::VectorXd &node1_ea, Eigen::VectorXd &node2_ea, Eigen::VectorXd &node3_ea, std::vector<double> &node_c_A, std::vector<double> &node_c_C, double dt_real, Eigen::VectorXd &Re, Eigen::VectorXd &Re_ea, Eigen::VectorXd &Re_Adif, Eigen::VectorXd &Re_Cdif, Eigen::Vector3d &MMe, int mat, Eigen::Vector3d &centroid, double time){

	// constants and parameters
	//double mu = 1.0e1; // [mili-Pa] material parameter
	double thetag_dot = 0.0; //0.0009; //0.005; // growth rate area_change/time_units, 0.01 is growth of 1% per second

	// example 1, biaxial test with no growth
	//double mu = 75000.0; // Pa, from a value we are using for skin 
	
	// MATERIAL CHOICE
	// stiffness 
	double mu = 0;
	if(mat==1){
		// cell stiffness 
		mu = 1000; // pN/ um^2 = Pa
	}else if(mat==2){
		// substrate is assummed much stiffer than the cell cytoskeleton
		mu = 1000000; // pN/um^2 = Pa
		//mu=1000;
		// if(time<5){
		// 	mu = 100000;
		// }else if(time<10){
		// 	mu = 20000;
		// }else if(time<15){
		// 	mu = 10000;
		// }else if(time<20){
		// 	mu = 5000;
		// }else{
		// 	mu = 1000;
		// }
	}
	// active stress
	double t_act = 0; 
	if(mat==1){
		 if(time<2){
		 	t_act = time*200.0; // Pa = pN/um^2 should be around 100-200 Pa on substrates < 100 kPa per: https://pubs.acs.org/doi/pdf/10.1021/am504135b 
		// }else if(time<5){
		// 	t_act = 200.0; 
		// }else if(time<6){
		// 	t_act = 200.0 + (time-5)*300.0;
		// }else if(time<10){
		// 	t_act = 500; 
		// }else if(time<11.0){
		// 	t_act = 500 + (time-10)*500;
		}else{
		 	t_act = 400;
		}
		// if(time<2){
		// 	t_act = time*200.0/2; // Pa = pN/um^2
		// }else{
		// 	t_act = 200;
		// }
	}else if(mat==2){
		t_act = 0;
	} 
	//double thetag_dot = 0.0; 
	//double thick = 0.005; // [meters] skin is ~0.5cm thick 
	double thick = 1; // um  

	// reference metric plus normal (reference meaning X_t, or configuration at time t)
	Eigen::Vector3d G1 = -1.0*node1_X + 1.0*node2_X;  // 3x1
	Eigen::Vector3d G2 = -1.0*node1_X + 1.0*node3_X; // 3x1
	Eigen::Vector3d G1xG2 = G1.cross(G2);
	Eigen::Vector3d N = G1xG2/G1xG2.norm(); // 3x1

	// interpolation of the euler almansi at the integration point 
	// this is ok for strains. 
	std::vector<Eigen::VectorXd> ea_t_nodes(3);
	ea_t_nodes[0]=node1_ea; 
	ea_t_nodes[1]=node2_ea;
	ea_t_nodes[2]=node3_ea; 
	Eigen::VectorXd ea_t_vec = 1.0/3.0*node1_ea + 1.0/3.0*node2_ea + 1.0/3.0*node3_ea;
	// the euler almansi strain is stored as vector. Need to reconstruct in Matrix
	Eigen::Matrix3d ea_t; 
	ea_t<<ea_t_vec(0),ea_t_vec(3),ea_t_vec(4),
		  ea_t_vec(3),ea_t_vec(1),ea_t_vec(5),
		  ea_t_vec(4),ea_t_vec(5),ea_t_vec(2);

	// will need to go from ea_t to be_t to get the area change and evaluate stress
	Eigen::Matrix3d Id; Id<<1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0;
	Eigen::Matrix3d aux = Id-2.0*ea_t;
	Eigen::Matrix3d be_t = aux.inverse();
	// note that the euler almansi strains that are stored correspond to the surface only
	// Moreover, they are stored with respect to the surface at the previous step 't'
	// which now has normal N 
	// the identity matrix used is the full 3D identity so be_t will have the extra 
	// normal component with no stretch. To get the surface only component we have. 
	Eigen::Matrix3d NdyadN = N*N.transpose();
	Eigen::Matrix3d be_t_s = be_t - NdyadN ;
	// do i need to store the thickness strain? I don't think so because the new
	// thickness strain is calculated based only on the current elastic deformation
	// we should just remember that the real be_t has some thickness strain to guarantee
	// incompressibility, just like the new be has to be adjusted for incompressibility 
	
	// deformed metric and normal (meaning the predicted state at time t+dt)
	Eigen::Vector3d g1 = -1.0*node1_x + 1.0*node2_x; // 3x1
	Eigen::Vector3d g2 = -1.0*node1_x + 1.0*node3_x; // 3x1
	Eigen::Vector3d g1crossg2 = g1.cross(g2);
	Eigen::Vector3d n = g1crossg2/g1crossg2.norm(); // 3x1
	Eigen::Matrix3d ndyadn = n*n.transpose();

	// put Gi in matrix G and get inverse
	Eigen::Matrix3d G; 
	G << G1(0),G2(0),N(0),
		 G1(1),G2(1),N(1),
		 G1(2),G2(2),N(2);
	Eigen::Matrix3d g; 
	g << g1(0),g2(0),n(0),
		 g1(1),g2(1),n(1),
		 g1(2),g2(2),n(2);
	//std::cout<<"G = "<<G<<"\ng="<<g<<"\n";
	Eigen::Matrix3d Ginv = G.inverse();
	Eigen::Vector3d G1inv; G1inv<<Ginv(0,0),Ginv(0,1),Ginv(0,2);
	Eigen::Vector3d G2inv; G2inv<<Ginv(1,0),Ginv(1,1),Ginv(1,2);
	// full deformation gradient at IP before static condensation 
	// this deformation gradient is the one corresponding to the time step, F_delta
	Eigen::Matrix3d F = g1*G1inv.transpose() + g2*G2inv.transpose()+n*N.transpose();
	//std::cout<<"F="<<F<<"\n";
	// not sure I need the following, shouldn't I be able to recast in terms of b?  
	Eigen::Matrix3d C = F.transpose()*F;
	Eigen::Matrix3d Cinv = C.inverse();
	Eigen::Matrix3d Cinv_s = Cinv - NdyadN;

	// gradient 
	Eigen::Matrix3d ginv = g.inverse();
	Eigen::Vector3d g1inv; g1inv<<ginv(0,0),ginv(0,1),ginv(0,2);
	Eigen::Vector3d g2inv; g2inv<<ginv(1,0),ginv(1,1),ginv(1,2);
	// gradA_xi = c_A[0]*[1,0] + c_A[1]*[0,1] + c_A[2]*[-1,-1]
	// gradA = ginvT*gradA_xi = c_A[0]*g1inv + c_A[1]*g2inv + c_A[2]*(-g1inv-g2inv)
	Eigen::Vector3d gradA = node_c_A[0]*g1inv + node_c_A[1]*g2inv + node_c_A[2]*(-g1inv-g2inv);
	Eigen::Vector3d gradC = node_c_C[0]*g1inv + node_c_C[1]*g2inv + node_c_C[2]*(-g1inv-g2inv);
	double c_A = (1./3.)*(node_c_A[0] + node_c_A[1] + node_c_A[2]);
	double c_C = (1./3.)*(node_c_C[0] + node_c_C[1] + node_c_C[2]);

	// Area changes to get the current normal strain 
	double theta_e_t = sqrt(be_t.determinant()); // use the total be_t with 1*NdyadN
	double theta_g = 1.0+thetag_dot*dt_real; // the growth rate is a constant  

	double theta = F.determinant(); // total area change from the deformation of the current step
	// the area changes satisfy theta*theta_e_t = theta_e*theta_g
	double theta_e = theta*theta_e_t/theta_g;
	// what about the volume change? Je_t = Je = 1, so Jdelta = Jg = theta_g

	// ABT 08/02/2018. Growth in response to elastic stretch
	//double K_thetae = 1.5;
	//double p_thetag = 0.0; // [1/sec]
	//if(theta_e>1.1){
	//	theta_g = 1.0+p_thetag*(theta_e/(K_thetae+theta_e));
	//}

	double lamda_N = 1.0/theta_e; // normal strain from the current Fe needed for incompressible
	double p = -mu*(lamda_N*lamda_N); // pressure for sigma_n = 0
	
	// STRESS refered to the 'current' configuration t, in my notes this is 'Sstar'
	Eigen::Matrix3d Spas = ((mu/theta_g)*be_t_s + theta_g*p*Cinv_s);
	// STRESS from myosin, i.e. active stress
	// isotropic for now 
	//Eigen::Matrix3d Sact = t_act*Id; 
	// fully anisotropic
	//Eigen::Matrix3d Sact = t_act*centroid*centroid.transpose(); 
	// somewhere in between isotropic (0.5) and anisootropic (0) 
	double kappa = 0.5;
	// need to re-scale because centroid is so close to zero, need unit vector from current
	// node to the centroid
	Eigen::Vector3d centroid_ei = 1./3.*(node1_x+node2_x+node3_x);
	Eigen::Vector3d a0 = centroid_ei-centroid;
	a0 = a0/a0.norm();
	Eigen::Matrix3d Sact = t_act*(kappa*Id +(1-2*kappa)*a0*a0.transpose()); 
	Eigen::Matrix3d S = Sact + Spas; 
	//std::cout<<"S_norm "<<S.norm()<<"\n";
	// UPDATE the strain 
	// well, the strain doesn't get updated right now, here we get the IP value of the updated
	// strain to interpolate it back to the nodes 
	Eigen::Matrix3d be = (1.0/(theta_g*theta_g))*F*be_t*F.transpose();
	Eigen::Matrix3d ea_s = 0.5*(Id-be.inverse());
	Eigen::VectorXd ea_s_vec(6);
	ea_s_vec<<ea_s(0,0),ea_s(1,1),ea_s(2,2),ea_s(0,1),ea_s(0,2),ea_s(1,2);

	// Given the stress fill out the element residual, lets leave the tangent alone
	Re.setZero();
	Re_ea.setZero();
	MMe.setZero();
	for(int i=0;i<3;i++){
	    for(int p=0;p<3;p++){
	        Eigen::Matrix3d temp; temp.setZero();
	        temp(p,i) = 1.0;
	        Eigen::Vector3d varg1 = -1.0*temp.col(0) + temp.col(1);
	        Eigen::Vector3d varg2 = -1.0*temp.col(0) + temp.col(2);
	        // varn =  cross(varg1,varg2)/norm(cross(varg1,varg2)) ??
	        Eigen::Matrix3d varF = varg1*G1inv.transpose() + varg2*G2inv.transpose(); //+outer(varn,Ginv[2,:])?
	        Eigen::Matrix3d varC = varF.transpose()*F + F.transpose()*varF;
	        Eigen::Matrix3d varE = 0.5*varC;
	        double trvarE = varE(0,0)+varE(1,1)+varE(2,2);
	        Eigen::VectorXd varE_voigt(6);
	        varE_voigt<<varE(0,0),varE(1,1),varE(2,2),2*varE(0,1),2*varE(0,2),2*varE(1,2);
	        Eigen::VectorXd S_voigt(6);
	        S_voigt<< S(0,0),S(1,1),S(2,2),S(0,1),S(0,2),S(1,2);
	        Re(i*3+p) = S_voigt.dot(varE_voigt)*G1xG2.norm()/2.0*thick;
	    }
	    for(int q=0;q<6;q++){
	    	//Re_ea(i*6+q) = (1./3.)*ea_s_vec(q)*G1xG2.norm()/2.0*thick;
	    	Re_ea(i*6+q) = (1./3.)*(ea_s_vec(q)-ea_t_nodes[i](q))*G1xG2.norm()/2.0*thick;
	    }
	    MMe(i) = (1./3.)*G1xG2.norm()/2.0*thick;
	    // diffusion 
	    Eigen::Vector3d  varN; varN.setZero();
	    varN(i) = 1.0;
	    Eigen::Vector3d gradNi = varN(0)*g1inv + varN(1)*g2inv + varN(2)*(-g1inv-g2inv);
	    // diffusion coefficient 
	    double D_A = 1.; // um^2/s fast diffusion of free integrins, can be higer
	    double D_C = 1e-10; 
	    // residual for the diffusion 
	    // int( (Cnew-Cold)/dt + gradC*gradN ) -> residual is 
	    // int( Cold + gradC*gradN )
	    Re_Adif(i)  = (c_A*(1./3.) - dt_real*D_A*gradA.dot(gradNi))*G1xG2.norm()/2.0*thick;
	    Re_Cdif(i)  = (c_C*(1./3.) - dt_real*D_C*gradC.dot(gradNi))*G1xG2.norm()/2.0*thick;
	}
}


void calculateRes(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Eigen::VectorXd> &ea_t,const NonDestructiveTriMesh &mesh, double dt_real, Eigen::VectorXd &Res, Eigen::VectorXd &Res_ea, Eigen::VectorXd &MM, int mat, double time){
	Res.setZero();
	Res_ea.setZero();
	MM.setZero();
	// START LOOP OVER ELEMENTS
	Eigen::VectorXd Re(9);
	Eigen::VectorXd Re_ea(18);
	Eigen::Vector3d MMe;
    int n_elem = mesh.num_triangles();
    for(int ei=0;ei<n_elem;ei++){
       	// connectivity of the linear elements
        int node1index = mesh.get_triangle(ei)[0];
    	int node2index = mesh.get_triangle(ei)[1];
    	int node3index = mesh.get_triangle(ei)[2];
        
        // previous time step nodal positions for this element
		Eigen::Vector3d node1_X; node1_X<<X_t[node1index][0],X_t[node1index][1],X_t[node1index][2];
		Eigen::Vector3d node2_X; node2_X<<X_t[node2index][0],X_t[node2index][1],X_t[node2index][2];
		Eigen::Vector3d node3_X; node3_X<<X_t[node3index][0],X_t[node3index][1],X_t[node3index][2];

		// previous time step euler almansi strains at the nodes
		Eigen::VectorXd node1_ea(6); 
		node1_ea<<ea_t[node1index](0),ea_t[node1index](1),ea_t[node1index](2),ea_t[node1index](3),ea_t[node1index](4),ea_t[node1index](5);
		Eigen::VectorXd node2_ea(6); 
		node2_ea<<ea_t[node2index](0),ea_t[node2index](1),ea_t[node2index](2),ea_t[node2index](3),ea_t[node2index](4),ea_t[node2index](5);
		Eigen::VectorXd node3_ea(6); 
		node3_ea<<ea_t[node3index](0),ea_t[node3index](1),ea_t[node3index](2),ea_t[node3index](3),ea_t[node3index](4),ea_t[node3index](5);

		// current nodal positions for this element
		Eigen::Vector3d node1_x; node1_x<<x_t[node1index][0],x_t[node1index][1],x_t[node1index][2];
		Eigen::Vector3d node2_x; node2_x<<x_t[node2index][0],x_t[node2index][1],x_t[node2index][2];
		Eigen::Vector3d node3_x; node3_x<<x_t[node3index][0],x_t[node3index][1],x_t[node3index][2];
						
        // and calculate the element Re
        Re.setZero();
        Re_ea.setZero();
        MMe.setZero();
				
        // subroutine to evaluate the element
        //std::cout<<"evaluating element "<<ei<<"\n";
        evalElementRe(node1_X,node2_X,node3_X,node1_x,node2_x,node3_x,node1_ea,node2_ea,node3_ea, dt_real, Re, Re_ea, MMe, mat, time);
        //std::cout<<"element residual norm "<<Re.norm()<<"\n";
        // LOOP OVER NODES
		for(int nodei=0;nodei<3;nodei++){
			// ASSEMBLE DISPLACEMENT RESIDUAL
			for(int coordi=0;coordi<3;coordi++){
				// residual
				Res(mesh.get_triangle(ei)[nodei]*3+coordi) += Re(nodei*3+coordi);
			}
			for(int coordi=0;coordi<6;coordi++){
				Res_ea(mesh.get_triangle(ei)[nodei]*6+coordi) += Re_ea(nodei*6+coordi);
			}
			MM(mesh.get_triangle(ei)[nodei]) += MMe(nodei);
		}
	}
}


void evalElementRe(Eigen::Vector3d &node1_X, Eigen::Vector3d &node2_X, Eigen::Vector3d &node3_X, Eigen::Vector3d &node1_x, Eigen::Vector3d &node2_x, Eigen::Vector3d &node3_x, Eigen::VectorXd &node1_ea, Eigen::VectorXd &node2_ea, Eigen::VectorXd &node3_ea, double dt_real, Eigen::VectorXd &Re, Eigen::VectorXd &Re_ea, Eigen::Vector3d &MMe, int mat, double time){

	// constants and parameters
	//double mu = 1.0e1; // [mili-Pa] material parameter
	double thetag_dot = 0.0; //0.0009; //0.005; // growth rate area_change/time_units, 0.01 is growth of 1% per second

	// example 1, biaxial test with no growth
	//double mu = 75000.0; // Pa, from a value we are using for skin 
	
	// MATERIAL CHOICE
	// stiffness 
	double mu = 0;
	if(mat==1){
		// cell stiffness 
		mu = 1000; // pN/ um^2 = Pa
	}else if(mat==2){
		// substrate is assummed much stiffer than the cell cytoskeleton
		mu = 1000000; // pN/um^2 = Pa
		//mu=1000;
		// if(time<5){
		// 	mu = 100000;
		// }else if(time<10){
		// 	mu = 20000;
		// }else if(time<15){
		// 	mu = 10000;
		// }else if(time<20){
		// 	mu = 5000;
		// }else{
		// 	mu = 1000;
		// }
	}
	// active stress
	double t_act = 0; 
	if(mat==1){
		 if(time<2){
		 	t_act = time*200.0; // Pa = pN/um^2 should be around 100-200 Pa on substrates < 100 kPa per: https://pubs.acs.org/doi/pdf/10.1021/am504135b 
		// }else if(time<5){
		// 	t_act = 200.0; 
		// }else if(time<6){
		// 	t_act = 200.0 + (time-5)*300.0;
		// }else if(time<10){
		// 	t_act = 500; 
		// }else if(time<11.0){
		// 	t_act = 500 + (time-10)*500;
		}else{
		 	t_act = 400;
		}
		// if(time<2){
		// 	t_act = time*200.0/2; // Pa = pN/um^2
		// }else{
		// 	t_act = 200;
		// }
	}else if(mat==2){
		t_act = 0;
	} 
	//double thetag_dot = 0.0; 
	//double thick = 0.005; // [meters] skin is ~0.5cm thick 
	double thick = 1; // um  

	// reference metric plus normal (reference meaning X_t, or configuration at time t)
	Eigen::Vector3d G1 = -1.0*node1_X + 1.0*node2_X;  // 3x1
	Eigen::Vector3d G2 = -1.0*node1_X + 1.0*node3_X; // 3x1
	Eigen::Vector3d G1xG2 = G1.cross(G2);
	Eigen::Vector3d N = G1xG2/G1xG2.norm(); // 3x1

	// interpolation of the euler almansi at the integration point 
	// this is ok for strains. 
	Eigen::VectorXd ea_t_vec = 1.0/3.0*node1_ea + 1.0/3.0*node2_ea + 1.0/3.0*node3_ea;
	//
	std::vector<Eigen::VectorXd> ea_t_nodes(3); 
	ea_t_nodes[0]=node1_ea;
	ea_t_nodes[1]=node2_ea; 
	ea_t_nodes[2]=node3_ea; 

	// the euler almansi strain is stored as vector. Need to reconstruct in Matrix
	Eigen::Matrix3d ea_t; 
	ea_t<<ea_t_vec(0),ea_t_vec(3),ea_t_vec(4),
		  ea_t_vec(3),ea_t_vec(1),ea_t_vec(5),
		  ea_t_vec(4),ea_t_vec(5),ea_t_vec(2);

	// will need to go from ea_t to be_t to get the area change and evaluate stress
	Eigen::Matrix3d Id; Id<<1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0;
	Eigen::Matrix3d aux = Id-2.0*ea_t;
	Eigen::Matrix3d be_t = aux.inverse();
	// note that the euler almansi strains that are stored correspond to the surface only
	// Moreover, they are stored with respect to the surface at the previous step 't'
	// which now has normal N 
	// the identity matrix used is the full 3D identity so be_t will have the extra 
	// normal component with no stretch. To get the surface only component we have. 
	Eigen::Matrix3d NdyadN = N*N.transpose();
	Eigen::Matrix3d be_t_s = be_t - NdyadN ;
	// do i need to store the thickness strain? I don't think so because the new
	// thickness strain is calculated based only on the current elastic deformation
	// we should just remember that the real be_t has some thickness strain to guarantee
	// incompressibility, just like the new be has to be adjusted for incompressibility 
	
	// deformed metric and normal (meaning the predicted state at time t+dt)
	Eigen::Vector3d g1 = -1.0*node1_x + 1.0*node2_x; // 3x1
	Eigen::Vector3d g2 = -1.0*node1_x + 1.0*node3_x; // 3x1
	Eigen::Vector3d g1crossg2 = g1.cross(g2);
	Eigen::Vector3d n = g1crossg2/g1crossg2.norm(); // 3x1
	Eigen::Matrix3d ndyadn = n*n.transpose();

	// put Gi in matrix G and get inverse
	Eigen::Matrix3d G; 
	G << G1(0),G2(0),N(0),
		 G1(1),G2(1),N(1),
		 G1(2),G2(2),N(2);
	Eigen::Matrix3d g; 
	g << g1(0),g2(0),n(0),
		 g1(1),g2(1),n(1),
		 g1(2),g2(2),n(2);
	//std::cout<<"G = "<<G<<"\ng="<<g<<"\n";
	Eigen::Matrix3d Ginv = G.inverse();
	Eigen::Vector3d G1inv; G1inv<<Ginv(0,0),Ginv(0,1),Ginv(0,2);
	Eigen::Vector3d G2inv; G2inv<<Ginv(1,0),Ginv(1,1),Ginv(1,2);
	// full deformation gradient at IP before static condensation 
	// this deformation gradient is the one corresponding to the time step, F_delta
	Eigen::Matrix3d F = g1*G1inv.transpose() + g2*G2inv.transpose()+n*N.transpose();
	//std::cout<<"F="<<F<<"\n";
	// not sure I need the following, shouldn't I be able to recast in terms of b?  
	Eigen::Matrix3d C = F.transpose()*F;
	Eigen::Matrix3d Cinv = C.inverse();
	Eigen::Matrix3d Cinv_s = Cinv - NdyadN;

	// Area changes to get the current normal strain 
	double theta_e_t = sqrt(be_t.determinant()); // use the total be_t with 1*NdyadN
	double theta_g = 1.0+thetag_dot*dt_real; // the growth rate is a constant  

	double theta = F.determinant(); // total area change from the deformation of the current step
	// the area changes satisfy theta*theta_e_t = theta_e*theta_g
	double theta_e = theta*theta_e_t/theta_g;
	// what about the volume change? Je_t = Je = 1, so Jdelta = Jg = theta_g

	// ABT 08/02/2018. Growth in response to elastic stretch
	//double K_thetae = 1.5;
	//double p_thetag = 0.0; // [1/sec]
	//if(theta_e>1.1){
	//	theta_g = 1.0+p_thetag*(theta_e/(K_thetae+theta_e));
	//}

	double lamda_N = 1.0/theta_e; // normal strain from the current Fe needed for incompressible
	double p = -mu*(lamda_N*lamda_N); // pressure for sigma_n = 0
	
	// STRESS refered to the 'current' configuration t, in my notes this is 'Sstar'
	Eigen::Matrix3d Spas = ((mu/theta_g)*be_t_s + theta_g*p*Cinv_s);
	// STRESS from myosin, i.e. active stress
	// isotropic for now 
	Eigen::Matrix3d Sact = t_act*Id; 
	Eigen::Matrix3d S = Sact + Spas; 
	//std::cout<<"S_norm "<<S.norm()<<"\n";
	// UPDATE the strain 
	// well, the strain doesn't get updated right now, here we get the IP value of the updated
	// strain to interpolate it back to the nodes 
	Eigen::Matrix3d be = (1.0/(theta_g*theta_g))*F*be_t*F.transpose();
	Eigen::Matrix3d ea_s = 0.5*(Id-be.inverse());
	Eigen::VectorXd ea_s_vec(6);
	ea_s_vec<<ea_s(0,0),ea_s(1,1),ea_s(2,2),ea_s(0,1),ea_s(0,2),ea_s(1,2);

	// Given the stress fill out the element residual, lets leave the tangent alone
	Re.setZero();
	Re_ea.setZero();
	MMe.setZero();
	for(int i=0;i<3;i++){
	    for(int p=0;p<3;p++){
	        Eigen::Matrix3d temp; temp.setZero();
	        temp(p,i) = 1.0;
	        Eigen::Vector3d varg1 = -1.0*temp.col(0) + temp.col(1);
	        Eigen::Vector3d varg2 = -1.0*temp.col(0) + temp.col(2);
	        // varn =  cross(varg1,varg2)/norm(cross(varg1,varg2)) ??
	        Eigen::Matrix3d varF = varg1*G1inv.transpose() + varg2*G2inv.transpose(); //+outer(varn,Ginv[2,:])?
	        Eigen::Matrix3d varC = varF.transpose()*F + F.transpose()*varF;
	        Eigen::Matrix3d varE = 0.5*varC;
	        double trvarE = varE(0,0)+varE(1,1)+varE(2,2);
	        Eigen::VectorXd varE_voigt(6);
	        varE_voigt<<varE(0,0),varE(1,1),varE(2,2),2*varE(0,1),2*varE(0,2),2*varE(1,2);
	        Eigen::VectorXd S_voigt(6);
	        S_voigt<< S(0,0),S(1,1),S(2,2),S(0,1),S(0,2),S(1,2);
	        Re(i*3+p) = S_voigt.dot(varE_voigt)*G1xG2.norm()/2.0*thick;
	    }
	    for(int q=0;q<6;q++){
	    	//Re_ea(i*6+q) = (1./3.)*ea_s_vec(q)*G1xG2.norm()/2.0*thick;
	    	Re_ea(i*6+q) = (1./3.)*(ea_s_vec(q)-ea_t_nodes[i](q))*G1xG2.norm()/2.0*thick;
	    }
	    MMe(i) = (1./3.)*G1xG2.norm()/2.0*thick;
	}
}



void calculateConstraintForce(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &x_t){
	// Force per node to guarantee attachment to the sphere (in the zebrafish sims)
	double dist;
	double radius = 350.0; // [m]
	double penalty_dist = 4000.0; // 1000.0; 1000 works
	F_constraint.setZero();
	int n_node = x_t.size();
	for(int ni=0;ni<n_node;ni++){
		dist = sqrt(x_t[ni][0]*x_t[ni][0] + x_t[ni][1]*x_t[ni][1] + x_t[ni][2]*x_t[ni][2]);
		// force should be proportional to the distance to the sphere of radius r 
		// and directed along the vector connecting the center of the sphere to the current x
		Vec3d u_x = x_t[ni]/dist;
		F_constraint(ni*3+0) = penalty_dist*(radius-dist)*u_x[0];
		F_constraint(ni*3+1) = penalty_dist*(radius-dist)*u_x[1];
		F_constraint(ni*3+2) = penalty_dist*(radius-dist)*u_x[2];
	}
}

void calculateForce_FA(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const std::vector<double> &c_C, std::vector<Vec3d> &fa_u, const NonDestructiveTriMesh &mesh,double kint){
	// force at a node depends directly on the displacement and FA stiffness 
	//std::cout<<"hello\n";
	//double k_FA = 1000; // pN/um
	double k_FA = kint;
	// FROM ANDRE, piece wise linear fit points
	int fFA_npts = 5; // for 0.1 nm/ns pull rate
	// int fFA_npts = 8; // for 2 and 10 nm/ns pull rate
	
	// OLD DATA 10 nm/ns pull rate
	// Eigen::VectorXd pxFA(fFA_npts); pxFA.setZero();
	// pxFA<< -1.4500e-06,  2.0542e-03,  4.0774e-03,  6.866e-03, 9.7555e-03;
	// Eigen::VectorXd pyFA(fFA_npts); pyFA.setZero();
	// pyFA<< 33.71792, 268.264 , 502.6462, 573.2455,503.2933;
	
	// 0.1 nm/ns pull rate
	// Eigen::VectorXd pxFA(fFA_npts); pxFA.setZero();
	// pxFA<< 0.        , 0.77360187e-03, 2.4288278e-03 , 5.27348466e-03, 6.473125e-03;
	// Eigen::VectorXd pyFA(fFA_npts); pyFA.setZero();
	// pyFA<< 6.18578872,  31.20492495, 257.81318168, 113.04735465, 218.57252843;
	
	// 2 nm/ns pull rate
	// Eigen::VectorXd pxFA(fFA_npts); pxFA.setZero();
	// pxFA<< -0.417375  ,  2.61887113,  6.17608203,  7.04473068, 10.44144885, 13.19262925, 16.69652913, 19.765875;
	// Eigen::VectorXd pyFA(fFA_npts); pyFA.setZero();
	// pyFA<< -20.89503266, 345.06755152, 227.63653033, 359.23102306, 197.25583105, 554.05871406, 384.17197124, 234.01261865;
	
	// 10 nm/ns pull rate
	// Eigen::VectorXd pxFA(fFA_npts); pxFA.setZero();
	// pxFA<< -0.05206917,  3.65738544,  5.08876181,  6.77812596, 11.61124615, 13.92569218, 14.53402364, 19.73099292;
	// Eigen::VectorXd pyFA(fFA_npts); pyFA.setZero();
	// pyFA<< 12.62828127, 593.67358827, 614.26466886, 619.23422115, 530.92336541, 545.8471786 , 688.95405979, 648.72872177;
	
	//std::cout<<"pxFA\n"<<pxFA;
	//std::cout<<"pyFA\n"<<pyFA;
	double rho_i_max = 100; // max integrin density in #/um^2 from "Influence of type I collagen surface density on fibroblast spreading, motility, and contractility.
	F_constraint.setZero();
	int n_node = x_t.size();
	for(int ni=0;ni<n_node;ni++){
		// force should be proportional to the displacement and weighted by the c_C field which has information about 
		// whether or not the node has a FA 
		// need the area of the node 
		// Get the triangles inciden in this vertex
		double Asum = 0.0;
		std::vector<size_t> onering = mesh.m_vertex_to_triangle_map[ni];
		for(int ti=0;ti<onering.size();ti++){
			// get the normal and area of each tringle
        	const Vec3st& tri = mesh.get_triangle(onering[ti]);
        	Eigen::Vector3d tri_n0; tri_n0<<x_t[tri[0]][0],x_t[tri[0]][1],x_t[tri[0]][2] ;
        	Eigen::Vector3d tri_n1; tri_n1<<x_t[tri[1]][0],x_t[tri[1]][1],x_t[tri[1]][2] ;
        	Eigen::Vector3d tri_n2; tri_n2<<x_t[tri[2]][0],x_t[tri[2]][1],x_t[tri[2]][2] ;
        	Eigen::Vector3d E1 = tri_n1 - tri_n0; 
        	Eigen::Vector3d E2 = tri_n2 - tri_n0; 
			Eigen::Vector3d normal_ti = E1.cross(E2);
			double area_ti = 0.5*normal_ti.norm();
			Asum += area_ti/3; 
		}
		// ABT feb 5, 2023
		// before I was updating the total displacement based on the deformation of the cell
		// in this timestep, assuming substrate fixed. But now doing this update before calling
		// this function, so fa should have the correct displacement of integrin 
		//Vec3d u_x = x_t[ni]-X_t[ni]+fa_u[ni];
		Vec3d u_x = fa_u[ni];
		F_constraint(ni*3+0) = -c_C[ni]*k_FA*u_x[0]*Asum*rho_i_max;
		F_constraint(ni*3+1) = -c_C[ni]*k_FA*u_x[1]*Asum*rho_i_max;
		F_constraint(ni*3+2) = -c_C[ni]*k_FA*u_x[2]*Asum*rho_i_max;
		// with difussion now values in absolute so no need to multiply by rho_i_max
		// F_constraint(ni*3+0) = -c_C[ni]*k_FA*u_x[0]*Asum;
		// F_constraint(ni*3+1) = -c_C[ni]*k_FA*u_x[1]*Asum;
		// F_constraint(ni*3+2) = -c_C[ni]*k_FA*u_x[2]*Asum;
		
		// with Andre data
		// Eigen::Vector3d Ffa; 
		// for(int fi=0;fi<fFA_npts-1;fi++){
			// for(int ci=0;ci<3;ci++){
				// if(abs(u_x[ci])>=pxFA(fi) && abs(u_x[ci])<pxFA(fi+1)){
					// double su = (u_x[ci]-pxFA(fi))/(pxFA(fi+1)-pxFA(fi));
					// if(u_x[ci]>0){
						// Ffa(ci) = (1-su)*pyFA(fi) + su*pxFA(fi+1);
					// }else{
						// Ffa(ci) = -(1-su)*pyFA(fi) - su*pxFA(fi+1);
					// }
				// }
			// }
		// }
		//std::cout<<"Ffa\n"<<Ffa;
		// F_constraint(ni*3+0) = -c_C[ni]*Ffa(0)*Asum*rho_i_max;
		// F_constraint(ni*3+1) = -c_C[ni]*Ffa(1)*Asum*rho_i_max;
		// F_constraint(ni*3+2) = -c_C[ni]*Ffa(2)*Asum*rho_i_max;

		// ABT feb 5, 2023, commenting out the update of fa vector
		// before I was updating the fa here from the cell movement, and then update the fa
		// from substrate movement in a different function, and dissipating energy in even 
		// one more function. Now trying to put all the fa updates into a single function
		// so comment out the update here...
		// update the reference position of the substrate, which doesn't move 
		// in the update FA I need to dissipate this 
		//fa_u[ni] = u_x;
	}
}

void calculateForce_FA_subs(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const NonDestructiveTriMesh &mesh, const std::vector<double> &c_C, const std::vector<Vec3d> &x_cell, std::vector<Vec3d> &fa_u, const NonDestructiveTriMesh &mesh_cell,double kint){
	F_constraint.setZero();
	// Loop over the nodes of the cell 
	int n_node = x_cell.size();
	for(int ni=0;ni<n_node;ni++){
		// get the location of the substrate point, not of the cell, from from cell to fa_u  
		Vec3d aux = x_cell[ni] - fa_u[ni];
		Eigen::Vector3d subs_X; subs_X<<aux[0],aux[1],aux[2];
		Eigen::Vector3d subs_x; subs_x<<aux[0],aux[1],aux[2];
		// find in which triangle this is, within the substrate mesh 
		// loop over the triangles of the substrate mesh 
		int n_elem = mesh.num_triangles();
	    for(int ei=0;ei<n_elem;ei++){
	       	// connectivity of the linear elements
	        int node1index = mesh.get_triangle(ei)[0];
	    	int node2index = mesh.get_triangle(ei)[1];
	    	int node3index = mesh.get_triangle(ei)[2];
	        
	        // previous time step nodal positions for this element
			Eigen::Vector3d node1_X; node1_X<<X_t[node1index][0],X_t[node1index][1],X_t[node1index][2];
			Eigen::Vector3d node2_X; node2_X<<X_t[node2index][0],X_t[node2index][1],X_t[node2index][2];
			Eigen::Vector3d node3_X; node3_X<<X_t[node3index][0],X_t[node3index][1],X_t[node3index][2];
			// new nodal positions
			Eigen::Vector3d node1_x; node1_x<<x_t[node1index][0],x_t[node1index][1],x_t[node1index][2];
			Eigen::Vector3d node2_x; node2_x<<x_t[node2index][0],x_t[node2index][1],x_t[node2index][2];
			Eigen::Vector3d node3_x; node3_x<<x_t[node3index][0],x_t[node3index][1],x_t[node3index][2];
			
			// check if point subs_s was in this triangle by getting barycentric coordinates
			//Eigen::Vector3d edge1 = node2_X - node1_X;
			//Eigen::Vector3d edge2 = node3_X - node1_X;
			Eigen::Vector3d edge1 = node2_x - node1_x;
			Eigen::Vector3d edge2 = node3_x - node1_x;
			Eigen::Vector3d normal = edge1.cross(edge2);
			double A_tri = 0.5*normal.norm();
			//Eigen::Vector3d edge_p = subs_X - node1_X; 
			Eigen::Vector3d edge_p = subs_x - node1_x; 
			Eigen::Vector3d A3n = edge1.cross(edge_p);
			double A3=-1;
			if(A3n(2)>=0){A3 = 0.5*A3n.norm();}
			double A2 = -1;
			Eigen::Vector3d A2n = edge_p.cross(edge2);
			if(A2n(2)>=0){A2 = 0.5*A2n.norm();}
			double A1 = A_tri - A2 - A3;
			if(A1>=0 && A2>=0 && A3>=0){
				//std::cout<<"point "<<subs_X(0)<<", "<<subs_X(1)<<"; is in tringle \n";
				//std::cout<<node1_X(0)<<", "<<node1_X(1)<<"; "<<node2_X(0)<<", "<<node2_X(1)<<"; "<<node3_X(0)<<", "<<node3_X(1)<<"\n"; 
				//std::cout<<"with barycentric coords "<<A1/A_tri<<", "<<A2/A_tri<<", "<<A3/A_tri<<"\n";
				// point is in triangle, we can actually see where this point has moved and get force
				// new x_subs position
				Vec3d subs_u = (A1/A_tri)*(x_t[node1index]-X_t[node1index])+(A2/A_tri)*(x_t[node2index]-X_t[node2index])+(A3/A_tri)*(x_t[node3index]-X_t[node3index]);
				// ABT Feb 2023. correct, assume that fa_u has the correct relative displacement between integrin and subs
				// calculate force with old displacement of the FA to keep the force balance between the two
				Vec3d u_x = fa_u[ni];
				// to calculate force we actually need node area 
				double Asum = 0.0;
				std::vector<size_t> onering = mesh_cell.m_vertex_to_triangle_map[ni];
				for(int ti=0;ti<onering.size();ti++){
					// get the normal and area of each tringle
		        	const Vec3st& tri = mesh_cell.get_triangle(onering[ti]);
		        	Eigen::Vector3d tri_n0; tri_n0<<x_cell[tri[0]][0],x_cell[tri[0]][1],x_cell[tri[0]][2] ;
		        	Eigen::Vector3d tri_n1; tri_n1<<x_cell[tri[1]][0],x_cell[tri[1]][1],x_cell[tri[1]][2] ;
		        	Eigen::Vector3d tri_n2; tri_n2<<x_cell[tri[2]][0],x_cell[tri[2]][1],x_cell[tri[2]][2] ;
		        	Eigen::Vector3d E1 = tri_n1 - tri_n0; 
		        	Eigen::Vector3d E2 = tri_n2 - tri_n0; 
					Eigen::Vector3d normal_ti = E1.cross(E2);
					double area_ti = 0.5*normal_ti.norm();
					Asum += area_ti/3; 
				} 
				// calculate force and distribute according to barycentric coordinates 
				// double k_FA = 1000; // pN/um, 1000 pN/um = 1 pN/nm
				double k_FA = kint;
				// ANDRE data
				int fFA_npts = 5; // for 0.1 nm/ns pull rate
				// int fFA_npts = 8; // for 2 and 10 nm/ns pull rate
				
				// OLD DATA 10 nm/ns pull rate
				// Eigen::VectorXd pxFA(fFA_npts); pxFA.setZero();
				// pxFA<< -1.4500e-06,  2.0542e-03,  4.0774e-03,  6.866e-03, 9.7555e-03;
				// Eigen::VectorXd pyFA(fFA_npts); pyFA.setZero();
				// pyFA<< 33.71792, 268.264 , 502.6462, 573.2455,503.2933;
				
				// 0.1 nm/ns pull rate
				// Eigen::VectorXd pxFA(fFA_npts); pxFA.setZero();
				// pxFA<< 0.        , 0.77360187e-03, 2.4288278e-03 , 5.27348466e-03, 6.473125e-03;
				// Eigen::VectorXd pyFA(fFA_npts); pyFA.setZero();
				// pyFA<< 6.18578872,  31.20492495, 257.81318168, 113.04735465, 218.57252843;
				
				// 2 nm/ns pull rate
				// Eigen::VectorXd pxFA(fFA_npts); pxFA.setZero();
				// pxFA<< -0.417375e-03  ,  2.61887113e-03,  6.17608203e-03,  7.04473068e-03, 10.44144885e-03, 13.19262925e-03, 16.69652913e-03, 19.765875e-03;
				// Eigen::VectorXd pyFA(fFA_npts); pyFA.setZero();
				// pyFA<< -20.89503266, 345.06755152, 227.63653033, 359.23102306, 197.25583105, 554.05871406, 384.17197124, 234.01261865;
				
				// 10 nm/ns pull rate
				// Eigen::VectorXd pxFA(fFA_npts); pxFA.setZero();
				// pxFA<< -0.05206917e-03,  3.65738544e-03,  5.08876181e-03,  6.77812596e-03, 11.61124615e-03, 13.92569218e-03, 14.53402364e-03, 19.73099292e-03;
				// Eigen::VectorXd pyFA(fFA_npts); pyFA.setZero();
				// pyFA<< 12.62828127, 593.67358827, 614.26466886, 619.23422115, 530.92336541, 545.8471786 , 688.95405979, 648.72872177;
				
				//
				double rho_i_max = 100; 
				// Comment loop below for without MD data
				// Eigen::Vector3d Ffa; 
				// for(int fi=0;fi<fFA_npts-1;fi++){
					// for(int ci=0;ci<3;ci++){
						// if(abs(u_x[ci])>=pxFA(fi) && abs(u_x[ci])<pxFA(fi+1)){
							// double su = (u_x[ci]-pxFA(fi))/(pxFA(fi+1)-pxFA(fi));
							// if(u_x[ci]>0){
								// Ffa(ci) = (1-su)*pyFA(fi) + su*pxFA(fi+1);
							// }
							// else{
								// Ffa(ci) = -(1-su)*pyFA(fi) - su*pxFA(fi+1);	
							// }
						// }
					// }
				// }
				// without ANDRE data
				F_constraint(node1index*3+0) += c_C[ni]*k_FA*u_x[0]*(A1/A_tri)*Asum*rho_i_max;
				F_constraint(node1index*3+1) += c_C[ni]*k_FA*u_x[1]*(A1/A_tri)*Asum*rho_i_max;
				F_constraint(node1index*3+2) += c_C[ni]*k_FA*u_x[2]*(A1/A_tri)*Asum*rho_i_max;
				F_constraint(node2index*3+0) += c_C[ni]*k_FA*u_x[0]*(A2/A_tri)*Asum*rho_i_max;
				F_constraint(node2index*3+1) += c_C[ni]*k_FA*u_x[1]*(A2/A_tri)*Asum*rho_i_max;
				F_constraint(node2index*3+2) += c_C[ni]*k_FA*u_x[2]*(A2/A_tri)*Asum*rho_i_max;
				F_constraint(node3index*3+0) += c_C[ni]*k_FA*u_x[0]*(A3/A_tri)*Asum*rho_i_max;
				F_constraint(node3index*3+1) += c_C[ni]*k_FA*u_x[1]*(A3/A_tri)*Asum*rho_i_max;
				F_constraint(node3index*3+2) += c_C[ni]*k_FA*u_x[2]*(A3/A_tri)*Asum*rho_i_max;
				
				// with 2-step diffusion no need for rho_i_max because absolute
				// F_constraint(node1index*3+0) += c_C[ni]*k_FA*u_x[0]*(A1/A_tri)*Asum;
				// F_constraint(node1index*3+1) += c_C[ni]*k_FA*u_x[1]*(A1/A_tri)*Asum;
				// F_constraint(node1index*3+2) += c_C[ni]*k_FA*u_x[2]*(A1/A_tri)*Asum;
				// F_constraint(node2index*3+0) += c_C[ni]*k_FA*u_x[0]*(A2/A_tri)*Asum;
				// F_constraint(node2index*3+1) += c_C[ni]*k_FA*u_x[1]*(A2/A_tri)*Asum;
				// F_constraint(node2index*3+2) += c_C[ni]*k_FA*u_x[2]*(A2/A_tri)*Asum;
				// F_constraint(node3index*3+0) += c_C[ni]*k_FA*u_x[0]*(A3/A_tri)*Asum;
				// F_constraint(node3index*3+1) += c_C[ni]*k_FA*u_x[1]*(A3/A_tri)*Asum;
				// F_constraint(node3index*3+2) += c_C[ni]*k_FA*u_x[2]*(A3/A_tri)*Asum;

				// with ANDRE data
				// F_constraint(node1index*3+0) += c_C[ni]*Ffa(0)*(A1/A_tri)*Asum*rho_i_max;
				// F_constraint(node1index*3+1) += c_C[ni]*Ffa(1)*(A1/A_tri)*Asum*rho_i_max;
				// F_constraint(node1index*3+2) += c_C[ni]*Ffa(2)*(A1/A_tri)*Asum*rho_i_max;
				// F_constraint(node2index*3+0) += c_C[ni]*Ffa(0)*(A2/A_tri)*Asum*rho_i_max;
				// F_constraint(node2index*3+1) += c_C[ni]*Ffa(1)*(A2/A_tri)*Asum*rho_i_max;
				// F_constraint(node2index*3+2) += c_C[ni]*Ffa(2)*(A2/A_tri)*Asum*rho_i_max;
				// F_constraint(node3index*3+0) += c_C[ni]*Ffa(0)*(A3/A_tri)*Asum*rho_i_max;
				// F_constraint(node3index*3+1) += c_C[ni]*Ffa(1)*(A3/A_tri)*Asum*rho_i_max;
				// F_constraint(node3index*3+2) += c_C[ni]*Ffa(2)*(A3/A_tri)*Asum*rho_i_max;

				// ABT feb 5, 2023
				// no more update of fa_u in this function, all updates of fa_u in the update_FA function
				// UPDATE the reference position wrt cell, which doesn't move in this update
				//u_x = fa_u[ni] - subs_u;
				//fa_u[ni] = u_x;
			}
		}
	}
}

void calculateForce_AP(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const std::vector<double> &c_A, const NonDestructiveTriMesh &mesh, double dt_real){
	// actin polymerization force should be applied on the boundary of the domain, should push outward
	// in the direction of the normal, with some random fluctuations 
	// initial random seed
	srand (time(NULL));
	// Loop over the nodes and find boundary 
	int n_node = X_t.size();
	for(int ni=0;ni<n_node;ni++){
		// check if this vertex is on the boundary
		if(mesh.m_is_boundary_vertex[ni]){
			// get the vector in Eigen format 
			Eigen::Vector3d ni_x; ni_x<<x_t[ni][0],x_t[ni][1],x_t[ni][2];
			// get the triangles incident on this vertex 
			std::vector<size_t> onering = mesh.m_vertex_to_triangle_map[ni];

			// loop over the one ring, there should be only two nodes in the one ring
			// on the boundary, so only two vectors expected
			int ct_aux = 0;
			Eigen::Vector3d ni_j1_x; 
			Eigen::Vector3d ni_j2_x; 
			int n1check = 0;
			int n2check = 0;
			// loop over one ring triangles
			for(int tj=0;tj<onering.size();tj++){
				// check if any of the vertex of the triangle and not ni
				// are on the boundary 
				const Vec3st& tri = mesh.get_triangle(onering[tj]);
				for(int enj=0;enj<3;enj++){
					int nj = tri[enj];
					if(mesh.m_is_boundary_vertex[nj] && nj!=ni){
						ct_aux += 1;
						// check if this should be the first or second node along 
						// the boundary curve, if cross product of boundary vector
						// and inside edge has positive 'z' then it is first vertex
						Eigen::Vector3d pb; pb<<x_t[nj][0],x_t[nj][1],x_t[nj][2];
						// actually need to know which is the interior node
						int nint;
						for(int enk=0;enk<3;enk++){
							if(tri[enk]!=nj && tri[enk]!=ni){nint=tri[enk];break;}
						}
						Eigen::Vector3d pint; pint<<x_t[nint][0],x_t[nint][1],x_t[nint][2];
						Eigen::Vector3d ncheck = (ni_x-pb).cross(pint-ni_x);
						if(ncheck[2]>0){
							ni_j1_x << x_t[nj][0],x_t[nj][1],x_t[nj][2];
							n1check = 1;
						}else{
							ni_j2_x << x_t[nj][0],x_t[nj][1],x_t[nj][2];
							n2check = 1;
						}
					}
				}
			}
			if(ct_aux>2){std::cout<<"boundary vertex with "<<ct_aux<<" neighbors\n";}
			if(n1check==0){std::cout<<"did not find the initial point of the boundary curve\n";}
			if(n2check==0){std::cout<<"did not find the final point of the boundary curve\n";}
			// given node ni and its neighbors on the boundary, get outward facing normal
			Eigen::Vector3d e3; e3<<0,0,1;
			Eigen::Vector3d normal_1 = -e3.cross(ni_x-ni_j1_x);
			normal_1 = normal_1/normal_1.norm();
			Eigen::Vector3d normal_2 = -e3.cross(ni_j2_x-ni_x);
			normal_2 = normal_2/normal_2.norm();
			double length1 = (ni_x-ni_j1_x).norm();
			double length2 = (ni_j2_x-ni_x).norm();
			double total_length = length1 + length2; 
			Eigen::Vector3d normal_i = length1/total_length*normal_1 + length2/total_length*normal_2;
			// add a random force in this direction 
			//int  factin_scale = (rand() % 500 )*total_length/2.0; // pN 
			// sample the factin as a poisson process 
			// https://en.wikipedia.org/wiki/Poisson_point_process
			double Pr =((double)rand()/((double)RAND_MAX+1.));
			double  factin_scale = 0;
			// https://www.pnas.org/doi/10.1073/pnas.0501435102
			double APrate = 10; // [1/s] but this needs also to be multiplied by density over length
			//double APrate_tot = 100*APrate*total_length/2; // where 50 is the density
			for(int ii=0;ii<100*total_length/2;ii++){
				if(Pr<APrate*dt_real*exp(-APrate*dt_real)){
					// APdensity*APforce*length, single
					// single APforce is order 1-10pN, 
					// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4810308/pdf/JCB_201512019.pdf
					// https://www.pnas.org/doi/pdf/10.1073/pnas.0607052104 
					// concentrations are 1-500uM, so in a thin shell near boundary
					// on the order of the elongation distance of 1-10nm 
					// ~5-100 monomers per um 
				 	factin_scale += 5; 
				}
			}
			F_constraint(ni*3+0) = factin_scale*normal_i(0);
			F_constraint(ni*3+1) = factin_scale*normal_i(1);
			F_constraint(ni*3+2) = factin_scale*normal_i(2);
		}
	}
}


void calculateBendingForce2D(Eigen::VectorXd &F_bending, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const NonDestructiveTriMesh &mesh){
	// actin polymerization force should be applied on the boundary of the domain, should push outward
	// in the direction of the normal, with some random fluctuations 
	// Loop over the nodes and find boundary 
	int n_node = X_t.size();
	for(int ni=0;ni<n_node;ni++){
		// check if this vertex is on the boundary
		if(mesh.m_is_boundary_vertex[ni]){
			// get the vector in Eigen format 
			Eigen::Vector3d ni_x; ni_x<<x_t[ni][0],x_t[ni][1],x_t[ni][2];
			// get the triangles incident on this vertex 
			std::vector<size_t> onering = mesh.m_vertex_to_triangle_map[ni];

			// loop over the one ring, there should be only two nodes in the one ring
			// on the boundary, so only two vectors expected
			int ct_aux = 0;
			Eigen::Vector3d ni_j1_x; 
			Eigen::Vector3d ni_j2_x; 
			int n1check = 0;
			int n2check = 0;
			// loop over one ring triangles
			for(int tj=0;tj<onering.size();tj++){
				// check if any of the vertex of the triangle and not ni
				// are on the boundary 
				const Vec3st& tri = mesh.get_triangle(onering[tj]);
				for(int enj=0;enj<3;enj++){
					int nj = tri[enj];
					if(mesh.m_is_boundary_vertex[nj] && nj!=ni){
						ct_aux += 1;
						// check if this should be the first or second node along 
						// the boundary curve, if cross product of boundary vector
						// and inside edge has positive 'z' then it is first vertex
						Eigen::Vector3d pb; pb<<x_t[nj][0],x_t[nj][1],x_t[nj][2];
						// actually need to know which is the interior node
						int nint;
						for(int enk=0;enk<3;enk++){
							if(tri[enk]!=nj && tri[enk]!=ni){nint=tri[enk];break;}
						}
						Eigen::Vector3d pint; pint<<x_t[nint][0],x_t[nint][1],x_t[nint][2];
						Eigen::Vector3d ncheck = (ni_x-pb).cross(pint-ni_x);
						if(ncheck[2]>0){
							ni_j1_x << x_t[nj][0],x_t[nj][1],x_t[nj][2];
							n1check = 1;
						}else{
							ni_j2_x << x_t[nj][0],x_t[nj][1],x_t[nj][2];
							n2check = 1;
						}
					}
				}
			}
			if(ct_aux>2){std::cout<<"boundary vertex with "<<ct_aux<<" neighbors\n";}
			if(n1check==0){std::cout<<"did not find the initial point of the boundary curve\n";}
			if(n2check==0){std::cout<<"did not find the final point of the boundary curve\n";}
			// get turning angle 
			Eigen::Vector3d ed1 = ni_x-ni_j1_x;
			ed1 = ed1/ed1.norm();
			Eigen::Vector3d ed2 = ni_j2_x-ni_x;
			ed2 = ed2/ed2.norm();
			double costurn = ed1.dot(ed2);
			Eigen::Vector3d e3sign = ed1.cross(ed2);
			double turna=0;
			if(e3sign(2)>0){
				 turna = costurn-1;
			}else{
				 turna = 1-costurn;
			}
			// calculate curvature
			double length1 = (ni_x-ni_j1_x).norm();
			double length2 = (ni_j2_x-ni_x).norm();
			double total_length = length1 + length2; 
			double kappa = turna*2/total_length;
			//std::cout<<"kappa "<<kappa<<"\n";

			// given node ni and its neighbors on the boundary, get outward facing normal
			Eigen::Vector3d e3; e3<<0,0,1;
			Eigen::Vector3d normal_1 = -e3.cross(ni_x-ni_j1_x);
			normal_1 = normal_1/normal_1.norm();
			Eigen::Vector3d normal_2 = -e3.cross(ni_j2_x-ni_x);
			normal_2 = normal_2/normal_2.norm();
			Eigen::Vector3d normal_i = length1/total_length*normal_1 + length2/total_length*normal_2;

			// force should be such that it minimizes curvature 
			double k_bend = 20.; // pN/um 
			F_bending(ni*3+0) = k_bend*kappa*normal_i(0);
			F_bending(ni*3+1) = k_bend*kappa*normal_i(1);
			F_bending(ni*3+2) = k_bend*kappa*normal_i(2);
			
		}
	}
}

Eigen::Vector4d calculateAreaTot(const std::vector<Vec3d> &x_t,const NonDestructiveTriMesh &mesh){
	int n_elem = mesh.num_triangles();
	double area_tot = 0;
	Eigen::Vector3d centroid; centroid.setZero();
    for(int ei=0;ei<n_elem;ei++){
		// connectivity of the linear elements
        int node1index = mesh.get_triangle(ei)[0];
    	int node2index = mesh.get_triangle(ei)[1];
    	int node3index = mesh.get_triangle(ei)[2];
        // nodal positions for this element
		Eigen::Vector3d node1_x; node1_x<<x_t[node1index][0],x_t[node1index][1],x_t[node1index][2];
		Eigen::Vector3d node2_x; node2_x<<x_t[node2index][0],x_t[node2index][1],x_t[node2index][2];
		Eigen::Vector3d node3_x; node3_x<<x_t[node3index][0],x_t[node3index][1],x_t[node3index][2];
		//
		Eigen::Vector3d area_ei_vec = 0.5*(node2_x-node1_x).cross(node3_x-node1_x);
		// centroid of this element 
		Eigen::Vector3d centroid_ei = 1./3.*(node1_x+node2_x+node3_x);

		double area_ei = area_ei_vec.norm();
		centroid += centroid_ei*area_ei;
		area_tot += area_ei;
    }
    centroid = centroid/area_tot;
    Eigen::Vector4d area; area<<area_tot,centroid(0),centroid(1),centroid(2);

    return area;
}


void calculateAreaConstraintForce(Eigen::VectorXd &F_areaConstraint, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const NonDestructiveTriMesh &mesh, double tot_area){
	int n_node = X_t.size();
	for(int ni=0;ni<n_node;ni++){
		// check if this vertex is on the boundary
		if(mesh.m_is_boundary_vertex[ni]){
			
			// get the vector in Eigen format 
			Eigen::Vector3d ni_x; ni_x<<x_t[ni][0],x_t[ni][1],x_t[ni][2];
			//std::cout<<"Node "<<ni<<"= "<<ni_x[0]<<", "<<ni_x[1]<<"; is on the boundary\n";
			// get the triangles incident on this vertex 
			std::vector<size_t> onering = mesh.m_vertex_to_triangle_map[ni];

			// loop over the one ring, there should be only two nodes in the one ring
			// on the boundary, so only two vectors expected
			int ct_aux = 0;
			Eigen::Vector3d ni_j1_x; 
			Eigen::Vector3d ni_j2_x; 
			int n1check = 0;
			int n2check = 0;
			// loop over one ring triangles
			for(int tj=0;tj<onering.size();tj++){
				// check if any of the vertex of the triangle and not ni
				// are on the boundary 
				const Vec3st& tri = mesh.get_triangle(onering[tj]);
				for(int enj=0;enj<3;enj++){
					int nj = tri[enj];
					if(mesh.m_is_boundary_vertex[nj] && nj!=ni){
						ct_aux += 1;
						//std::cout<<"found the #"<<ct_aux<<" neighbor on boundary, node "<<nj<<"\n";
						//std::cout<<"it is on triangle "<<tj<<"\n";
						// check if this should be the first or second node along 
						// the boundary curve, if cross product of boundary vector
						// and inside edge has positive 'z' then it is first vertex
						Eigen::Vector3d pb; pb<<x_t[nj][0],x_t[nj][1],x_t[nj][2];
						// actually need to know which is the interior node
						int nint;
						for(int enk=0;enk<3;enk++){
							if(tri[enk]!=nj && tri[enk]!=ni){nint=tri[enk];break;}
						}
						//std::cout<<"the third (interior) node of the triangle is "<<nint<<"\n";
						Eigen::Vector3d pint; pint<<x_t[nint][0],x_t[nint][1],x_t[nint][2];
						Eigen::Vector3d ncheck = (ni_x-pb).cross(pint-ni_x);
						if(ncheck[2]>0){
							ni_j1_x << x_t[nj][0],x_t[nj][1],x_t[nj][2];
							n1check = 1;
							//std::cout<<"neighbor #"<<ct_aux<<" is the beginning of boundary curve segment\n";
						}else{
							ni_j2_x << x_t[nj][0],x_t[nj][1],x_t[nj][2];
							n2check = 1;
							//std::cout<<"neighbor "<<ct_aux<<" is the end of boundary curve segment\n";
						}
					}
				}
			}
			if(ct_aux>2){std::cout<<"boundary vertex with "<<ct_aux<<" neighbors\n";}
			if(n1check==0){std::cout<<"did not find the initial point of the boundary curve\n";}
			if(n2check==0){std::cout<<"did not find the final point of the boundary curve\n";}
			// get turning angle 
			Eigen::Vector3d ed1 = ni_x-ni_j1_x;
			ed1 = ed1/ed1.norm();
			Eigen::Vector3d ed2 = ni_j2_x-ni_x;
			ed2 = ed2/ed2.norm();
			double costurn = ed1.dot(ed2);
			Eigen::Vector3d e3sign = ed1.cross(ed2);
			double turna=0;
			if(e3sign(2)>0){
				 turna = acos(costurn);
			}else{
				turna = -acos(costurn);
			}
			// calculate curvature
			double length1 = (ni_x-ni_j1_x).norm();
			double length2 = (ni_j2_x-ni_x).norm();
			double total_length = length1 + length2; 
			double kappa = (3.14159 - turna)*2/total_length;

			// given node ni and its neighbors on the boundary, get outward facing normal
			Eigen::Vector3d e3; e3<<0,0,1;
			Eigen::Vector3d normal_1 = -e3.cross(ni_x-ni_j1_x);
			normal_1 = normal_1/normal_1.norm();
			Eigen::Vector3d normal_2 = -e3.cross(ni_j2_x-ni_x);
			normal_2 = normal_2/normal_2.norm();
			Eigen::Vector3d normal_i = length1/total_length*normal_1 + length2/total_length*normal_2;

			// force should be such that it minimizes curvature 
			double k_area = 1; // pN/um 
			double area_0 = 10*10*3.1416;
			F_areaConstraint(ni*3+0) = -k_area*(tot_area-area_0)/2*total_length*normal_i(0);
			F_areaConstraint(ni*3+1) = -k_area*(tot_area-area_0)/2*total_length*normal_i(1);
			F_areaConstraint(ni*3+2) = -k_area*(tot_area-area_0)/2*total_length*normal_i(2);
			
		}
	}
}

void calculatePressureForce(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &x_t, const Eigen::VectorXd &MM,const NonDestructiveTriMesh &mesh, double t_tot ){
	F_constraint.setZero();
	double pres = 2000.0; // [Pa] 
	if(t_tot<1){
		pres=2000.0*t_tot;
		//pres = 0.0;
	}
	if(t_tot>2){
		pres=2000.0*(3-t_tot);
		//pres = 0.0;
	}
	if(t_tot>3){
		pres=0.0;
	}
	
	// loop over nodes and get the normal vector
	double thick = 0.1; //[m] this should be passed actually, I just know it is 0.1 for everywhere 
	for(int ni=0;ni<x_t.size();ni++){
		Eigen::Vector3d normal; normal.setZero();
		// Get the triangles inciden in this vertex
		double Asum = 0.0;
		std::vector<size_t> onering = mesh.m_vertex_to_triangle_map[ni];
		for(int ti=0;ti<onering.size();ti++){
			// get the normal and area of each tringle
        	const Vec3st& tri = mesh.get_triangle(onering[ti]);
        	Eigen::Vector3d tri_n0; tri_n0<<x_t[tri[0]][0],x_t[tri[0]][1],x_t[tri[0]][2] ;
        	Eigen::Vector3d tri_n1; tri_n1<<x_t[tri[1]][0],x_t[tri[1]][1],x_t[tri[1]][2] ;
        	Eigen::Vector3d tri_n2; tri_n2<<x_t[tri[2]][0],x_t[tri[2]][1],x_t[tri[2]][2] ;
        	Eigen::Vector3d E1 = tri_n1 - tri_n0; 
        	Eigen::Vector3d E2 = tri_n2 - tri_n0; 
			Eigen::Vector3d normal_ti = E1.cross(E2);
			double area_ti = 0.5*normal_ti.norm();
			normal_ti = normal_ti/normal_ti.norm();
			
			// add to the average calculation
			normal += area_ti*normal_ti;
			Asum += area_ti; 
		}
		normal = normal/Asum;
		normal = normal/normal.norm();
		
		//std::cout<<"vertex "<<ni<<", normal = "<<normal(0)<<", "<<normal(1)<<", "<<normal(2)<<", area = "<<Asum<<"\n";
		// pressure is constant traction, needs to be scaled by area of each node
		F_constraint(ni*3+0) = -pres*(MM(ni)/thick)*normal(0);
		F_constraint(ni*3+1) = -pres*(MM(ni)/thick)*normal(1);
		F_constraint(ni*3+2) = -pres*(MM(ni)/thick)*normal(2);
		//F_constraint(ni*3+0) = 0.0;
		//F_constraint(ni*3+1) = 0.0;
		//F_constraint(ni*3+2) = -1.0;
		//std::cout<<"F_constraint = "<<F_constraint(ni*3+0)<<", "<<F_constraint(ni*3+1)<<", "<<F_constraint(ni*3+2)<<"\n";  
	}
}

void calculateBendingRes(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const NonDestructiveTriMesh &mesh, Eigen::VectorXd &Res_bending){
	// Given the coordinates of the nodes at the end of time step t, and the new positions at t+dt, get the
	// forces for the next time step, and that includes some bending opposing force. 
	// this is mostly to regularize the deformation and avoid the bad crumpling effect

	Res_bending.setZero();
	// loop over edges
	int n_edges = mesh.m_edges.size();
	//std::cout<<"loop over "<<n_edges<<"edges\n";
	for(int e=0;e<n_edges;e++){
		// if it is an interior triangle then get the triangles incident on it
		if(!mesh.m_is_boundary_edge[e]){
			// get triangles incident to the edge
        	int t0 = mesh.m_edge_to_triangle_map[e][0];
        	int t1 = mesh.m_edge_to_triangle_map[e][1];
        	const Vec3st& tri0 = mesh.get_triangle(t0);
        	const Vec3st& tri1 = mesh.get_triangle(t1);
			// get vertex opposite the edge for each triangle
        	int v2_t0 = mesh.get_third_vertex( e, tri0 );
        	int v2_t1 = mesh.get_third_vertex( e, tri1 );
			// get also the two vertices that form this edge
			int v0 = mesh.m_edges[e][0];      
			int v1 = mesh.m_edges[e][1];

			// get the reference coordinates of the vertices
			Eigen::Vector3d v0_X, v1_X, v2t0_X, v2t1_X;
			v0_X<<X_t[v0][0],X_t[v0][1],X_t[v0][2];
			v1_X<<X_t[v1][0],X_t[v1][1],X_t[v1][2];
			v2t0_X<<X_t[v2_t0][0],X_t[v2_t0][1],X_t[v2_t0][2];
			v2t1_X<<X_t[v2_t1][0],X_t[v2_t1][1],X_t[v2_t1][2];
			// get the deformed coordinates of the vertices
			Eigen::Vector3d v0_x, v1_x, v2t0_x, v2t1_x;
			v0_x<<x_t[v0][0],x_t[v0][1],x_t[v0][2];
			v1_x<<x_t[v1][0],x_t[v1][1],x_t[v1][2];
			v2t0_x<<x_t[v2_t0][0],x_t[v2_t0][1],x_t[v2_t0][2];
			v2t1_x<<x_t[v2_t1][0],x_t[v2_t1][1],x_t[v2_t1][2];

			// put the indices and coordinates of the vertices in vectors
			std::vector<int> v_vec_ind(4,0);
			v_vec_ind[0]= v0; v_vec_ind[1]= v1;
			v_vec_ind[2]= v2_t0; v_vec_ind[3]= v2_t1;
			std::vector<Eigen::Vector3d> v_vec_X(4,Eigen::Vector3d());
			v_vec_X[0] = v0_X; v_vec_X[1] = v1_X; 
			v_vec_X[2] = v2t0_X; v_vec_X[3] = v2t1_X; 
			std::vector<Eigen::Vector3d> v_vec_x(4,Eigen::Vector3d());
			v_vec_x[0] = v0_x; v_vec_x[1] = v1_x; 
			v_vec_x[2] = v2t0_x; v_vec_x[3] = v2t1_x; 

			// ENERGY from bending
			double Energy_b = calculateBendingEnergy(v_vec_X,v_vec_x,false); 
			//std::cout<<"Energy_b = "<<Energy_b<<"\n";
			// variation of the bending in terms of the displacement to update force
			// loop over the four nodes that affect this edge
			double delta = 1e-6;
			double Eb_plus;
			for(int ni=0;ni<4;ni++){
				for(int ci=0;ci<3;ci++){
					// variation of vj, first numerically for simplicity
					v_vec_x[ni](ci) += delta; 
					Eb_plus = calculateBendingEnergy(v_vec_X,v_vec_x,false); 
					Res_bending(v_vec_ind[ni]*3+ci)+=(Eb_plus-Energy_b)/delta;
					v_vec_x[ni](ci) -= delta;
				}
			}
		}
	}

}

double calculateBendingEnergy(const std::vector<Eigen::Vector3d> &v_vec_X , const std::vector<Eigen::Vector3d> &v_vec_x, bool print_info){
	
	// beding energy is controlled by a material parameter
	//double mu_b = 10000.0;
	double mu_b = 10.0; //[Pa]
	double h = 0.1; // [m], thickness

	// unpack
	Eigen::Vector3d v0_X, v1_X, v2t0_X, v2t1_X;
	Eigen::Vector3d v0_x, v1_x, v2t0_x, v2t1_x;
	v0_X = v_vec_X[0];
	v1_X = v_vec_X[1];
	v2t0_X = v_vec_X[2];
	v2t1_X = v_vec_X[3];
	v0_x = v_vec_x[0];
	v1_x = v_vec_x[1];
	v2t0_x = v_vec_x[2];
	v2t1_x = v_vec_x[3];

	// get the angle between the two triangles before, based on X
	// normal of triangles
	Eigen::Vector3d edge_X = v1_X-v0_X;
	Eigen::Vector3d edget0_X = v2t0_X-v0_X;
	Eigen::Vector3d Nt0_X = edge_X.cross(edget0_X);
	// area
	double At0 = 0.5*Nt0_X.norm();
	// normal
	Nt0_X = Nt0_X/Nt0_X.norm();
	//
	Eigen::Vector3d edget1_X = v2t1_X-v0_X;
	Eigen::Vector3d Nt1_X = edget1_X.cross(edge_X);
	// area
	double At1 = 0.5*Nt1_X.norm();
	Nt1_X = Nt1_X/Nt1_X.norm();
	// get the angle directly from normals
	double N0dotN1 = Nt0_X.dot(Nt1_X)>=1 ? 0.9999999999 : Nt0_X.dot(Nt1_X);
	double thetae_X = acos(N0dotN1);
	if(print_info){
	std::cout<<"initial areas, At0="<<At0<<", At1="<<At1<<"\n";
	std::cout<<"reference normals, N0=("<<Nt0_X(0)<<","<<Nt0_X(1)<<","<<Nt0_X(2)<<"), N1=("<<Nt1_X(0)<<","<<Nt0_X(1)<<","<<Nt0_X(2)<<")\n";
	std::cout<<"dot product = "<<Nt0_X.dot(Nt1_X)<<"\n";
	std::cout<<"reference angle, theta_X = "<<thetae_X<<"\n";}
	// do the same for the deformed triangle
	Eigen::Vector3d edge_x = v1_x-v0_x;
	Eigen::Vector3d edget0_x = v2t0_x-v0_x;
	Eigen::Vector3d Nt0_x = edge_x.cross(edget0_x);
	Nt0_x = Nt0_x/Nt0_x.norm();
	//
	Eigen::Vector3d edget1_x = v2t1_x-v0_x;
	Eigen::Vector3d Nt1_x = edget1_x.cross(edge_x);
	Nt1_x = Nt1_x/Nt1_x.norm();
	// deformed angle
	double n0dotn1 = Nt0_x.dot(Nt1_x)>=1 ? 0.9999999999 : Nt0_x.dot(Nt1_x);
	double thetae_x = acos(n0dotn1);
	if(print_info){
	std::cout<<"deformed normals, n0=("<<Nt0_x(0)<<","<<Nt0_x(1)<<","<<Nt0_x(2)<<"), n1=("<<Nt1_x(0)<<","<<Nt0_x(1)<<","<<Nt0_x(2)<<")\n";
	std::cout<<"dot product = "<<Nt0_x.dot(Nt1_x)<<"\n";
	std::cout<<"deformed angle, theta_x = "<<thetae_x<<"\n";}
	// length of the edge
	double lbar_edge = edge_X.norm();
	// barycentric length
	double lstar_edge = (At0+At1)/(3*lbar_edge);

	// ENERGY from bending [N*m]
	double Energy_b = mu_b*h*h*h*(thetae_x*thetae_x)*lbar_edge/lstar_edge; 


	return Energy_b;
}

//ABT Feb 5, 2023
// change this function a lot so that it can be called outside of 'calculateForces'
// because it is too many changes, going to create a new function with different arguments
// this is called function overloading in C++
void updateFA(std::vector<double> &c_C, std::vector<Vec3d> &fa_u, double dt_real, const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, double kint){
	// Update focal adhesion concentration 
	// Bell model of adhesion, no diffusion, point-wise 
	// October 2022, update to change reference configuration based on disolution of bonds
	// and creation of new ones  
	int n_node = c_C.size();
	for(int ni=0;ni<n_node;ni++){
		// calculate single integrin force directly from stiffness and displacement 
		//double k_FA = 1000; // pN/um 
		double k_FA = kint;
		double force_FA = k_FA* sqrt(fa_u[ni][0]*fa_u[ni][0] + fa_u[ni][1]*fa_u[ni][1]);
		// update thinking of some well potential with the force
		double k_on = 1; // 1/s 
		double ka = 0.004; // 1/s
		double kb = 10; // 1/s
		double Fa = 15; //pN 
		double Fb = 15; // pN based on Schumacher et al. 2021
		double k_off = ka*exp(force_FA/Fa) + kb*exp(-force_FA/Fb); 
		//std::cout<<"force per integrin "<<force_FA<<"\n";
		double c_C_old = c_C[ni];
		c_C[ni] = c_C[ni] + dt_real*k_on*(1-c_C[ni]) - dt_real*k_off*c_C[ni]; 
		if(c_C[ni]<0){c_C[ni]=0;}

		// now update the reference state
		// this is the displacement of the old bonds 
		// it accumulates from before 
		//Vec3d u_x = x_t[ni]-X_t[ni]+fa_u[ni];
		// but as some of these bonds dissapear and new ones are formed
		// the actual reference should be a weighted sum of the new reference and the old
		// the new ones have zero displacement, so just need to drift the old ones 
		//fa_u[ni] = u_x*(c_C_old - dt_real*k_off*c_C_old)/c_C[ni];
		// Alternatively, just dissipate here 
		if(c_C_old - dt_real*k_off*c_C_old>0 && c_C[ni]>0){
			fa_u[ni] = fa_u[ni]*(c_C_old - dt_real*k_off*c_C_old)/c_C[ni];
		}else{
			fa_u[ni] = 0.0*fa_u[ni];
		}
	}
}

// ABT feb 5, 2023
// NEW update FA function that can be called outside of 'calculateForces'
// function should take previous and new nodal positions of cell and substrate and use to calculate the 
// new fa and the new c_C
void updateFA(std::vector<double> &c_A, std::vector<double> &c_C, std::vector<Vec3d> &fa_u, double dt_real, const std::vector<Vec3d> &X_t_cell, const std::vector<Vec3d> &x_t_cell,const NonDestructiveTriMesh &mesh_cell, const std::vector<Vec3d> &X_t_subs, const std::vector<Vec3d> &x_t_subs,const NonDestructiveTriMesh &mesh_subs, double kint){
	// Update focal adhesion concentration 
	// Bell model of adhesion, no diffusion, point-wise 
	// October 2022, update to change reference configuration based on disolution of bonds
	// and creation of new ones  
	int n_node = c_C.size();
	for(int ni=0;ni<n_node;ni++){
		// for every node of the cell, need the displacement from the cell side but also the 
		// interpolated displacement from the substrate 
		//
		// displacement from cell side
		Vec3d cell_u = x_t_cell[ni]-X_t_cell[ni];
		// displacement from subs side: need to first find which triangle the X_cell was on, then interpolate subs_u
		Vec3d subs_u; 
		// get the location of the substrate point, not of the cell, from from cell to fa_u  
		Vec3d aux = X_t_cell[ni] - fa_u[ni];
		Eigen::Vector3d subs_X; subs_X<<aux[0],aux[1],aux[2];
		//Eigen::Vector3d subs_x; subs_x<<aux[0],aux[1],aux[2];
		// find in which triangle this is, within the substrate mesh 
		// loop over the triangles of the substrate mesh 
		int n_elem = mesh_subs.num_triangles();
	    for(int ei=0;ei<n_elem;ei++){
	       	// connectivity of the linear elements
	        int node1index = mesh_subs.get_triangle(ei)[0];
	    	int node2index = mesh_subs.get_triangle(ei)[1];
	    	int node3index = mesh_subs.get_triangle(ei)[2];
	        
	        // previous time step nodal positions for this element
			Eigen::Vector3d node1_X_subs; node1_X_subs<<X_t_subs[node1index][0],X_t_subs[node1index][1],X_t_subs[node1index][2];
			Eigen::Vector3d node2_X_subs; node2_X_subs<<X_t_subs[node2index][0],X_t_subs[node2index][1],X_t_subs[node2index][2];
			Eigen::Vector3d node3_X_subs; node3_X_subs<<X_t_subs[node3index][0],X_t_subs[node3index][1],X_t_subs[node3index][2];
			// new nodal positions
			Eigen::Vector3d node1_x_subs; node1_x_subs<<x_t_subs[node1index][0],x_t_subs[node1index][1],x_t_subs[node1index][2];
			Eigen::Vector3d node2_x_subs; node2_x_subs<<x_t_subs[node2index][0],x_t_subs[node2index][1],x_t_subs[node2index][2];
			Eigen::Vector3d node3_x_subs; node3_x_subs<<x_t_subs[node3index][0],x_t_subs[node3index][1],x_t_subs[node3index][2];
			
			// check if point subs_s was in this triangle by getting barycentric coordinates
			Eigen::Vector3d edge1 = node2_X_subs - node1_X_subs;
			Eigen::Vector3d edge2 = node3_X_subs - node1_X_subs;
			//Eigen::Vector3d edge1 = node2_x_subs - node1_x_subs;
			//Eigen::Vector3d edge2 = node3_x_subs - node1_x_subs;
			Eigen::Vector3d normal = edge1.cross(edge2);
			double A_tri = 0.5*normal.norm();
			Eigen::Vector3d edge_p = subs_X - node1_X_subs; 
			//Eigen::Vector3d edge_p = subs_x - node1_x_subs; 
			Eigen::Vector3d A3n = edge1.cross(edge_p);
			double A3=-1;
			if(A3n(2)>=0){A3 = 0.5*A3n.norm();}
			double A2 = -1;
			Eigen::Vector3d A2n = edge_p.cross(edge2);
			if(A2n(2)>=0){A2 = 0.5*A2n.norm();}
			double A1 = A_tri - A2 - A3;
			if(A1>=0 && A2>=0 && A3>=0){
				//std::cout<<"point "<<subs_X(0)<<", "<<subs_X(1)<<"; is in tringle \n";
				//std::cout<<node1_X(0)<<", "<<node1_X(1)<<"; "<<node2_X(0)<<", "<<node2_X(1)<<"; "<<node3_X(0)<<", "<<node3_X(1)<<"\n"; 
				//std::cout<<"with barycentric coords "<<A1/A_tri<<", "<<A2/A_tri<<", "<<A3/A_tri<<"\n";
				// point is in triangle, we can actually see where this point has moved and get force
				// new x_subs position
				subs_u = (A1/A_tri)*(x_t_subs[node1index]-X_t_subs[node1index])+(A2/A_tri)*(x_t_subs[node2index]-X_t_subs[node2index])+(A3/A_tri)*(x_t_subs[node3index]-X_t_subs[node3index]);
			}
		}

		// COMPUTE INTEGRIN-LIGAND FORCE 
		// so now we have the correct subs displacement and also the correct cell displacement 
		// use total displacement to get the updated fa_u, remember fa_u is a vector going from 
		// the substrate point to the cell point. Right now fa_u has the correct displacement between subs and cell from previous time step
		// but now there are new cell and subs locations, use those to update fa_u 
		fa_u[ni] = fa_u[ni] + cell_u - subs_u;  
		// with this new integrin-ligand distance compute the force
		// calculate single integrin force directly from stiffness and displacement 
		//double k_FA = 1000; // pN/um 
		double k_FA = kint;
		double force_FA = k_FA* sqrt(fa_u[ni][0]*fa_u[ni][0] + fa_u[ni][1]*fa_u[ni][1]);
		// update thinking of some well potential with the force
		double k_on = 1; // 1/s 
		double ka = 0.004; // 1/s
		double kb = 10; // 1/s
		double Fa = 15; //pN 
		double Fb = 15; // pN based on Schumacher et al. 2021
		double k_off = ka*exp(force_FA/Fa) + kb*exp(-force_FA/Fb); 
		//std::cout<<"force per integrin "<<force_FA<<"\n";
		double c_C_old = c_C[ni];
		
		c_C[ni] = c_C[ni] + dt_real*k_on*(1-c_C[ni]) - dt_real*k_off*c_C[ni]; 

		// New way of doing it with the 2-step model 
		//double c_Cmax = 100; // integrin bonds/um^2
		//c_C[ni] = c_C[ni] + dt_real*k_on*c_A[ni]*(c_Cmax-c_C[ni]) - dt_real*k_off*c_C[ni];
		//c_A[ni] = c_A[ni] - dt_real*k_on*c_A[ni]*(c_Cmax-c_C_old) + dt_real*k_off*c_C_old;

		if(c_C[ni]<0){c_C[ni]=0;}
		if(c_A[ni]<0){c_A[ni]=0;}

		// now one more update the reference state
		// this is the displacement of the old bonds 
		// it accumulates from before 
		//Vec3d u_x = fa_u + cell_u - subs_u ;
		// but as some of these bonds dissapear and new ones are formed
		// the actual reference should be a weighted sum of the new reference and the old
		// the new ones have zero displacement, they should form stress-free so just need to drift the old ones 
		//fa_u[ni] = u_x*(c_C_old - dt_real*k_off*c_C_old)/c_C[ni];
		// Alternatively, just dissipate here 
		if(c_C_old - dt_real*k_off*c_C_old>0 && c_C[ni]>0){
			fa_u[ni] = fa_u[ni]*(c_C_old - dt_real*k_off*c_C_old)/c_C[ni];
		}else{
			fa_u[ni] = 0.0*fa_u[ni];
		}
	}

}

