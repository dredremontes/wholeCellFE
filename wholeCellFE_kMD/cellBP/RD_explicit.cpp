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

// Explicit time step for reaction-diffusion system 

// For now, simply update the FA concentration, not based on force? Random? 
void updateC_Explicit(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t, std::vector<Eigen::VectorXd> &ea_t, 
	std::vector<double> &c_B, std::vector<double> &c_C, std::vector<double> &c_BC, double dt, double dt_real, const NonDestructiveTriMesh &mesh, double t_tot){

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
	calculateRes(X_t,x_t,ea_t,mesh,dt_real, Res,Res_ea,MM); 
	//******************************************************//

	// Calculate the force constraint, forces needed to enforce things like:
	// * growing tissue wants to stick to the sphere
	Eigen::VectorXd F_FA(n_node*3);F_FA.setZero();
	Eigen::VectorXd F_AP(n_node*3);F_AP.setZero();
	
	//calculateConstraintForce(F_constraint,x_t);

	// Force for focal adhesions, need previous and new positions
	// as well as whether there is a FA, maybe will need additional parameters
	// or state variables for each adhesion, but first just one
	calculateForce_FA(F_FA,X_t,x_t,c_B);

	// Also need to add some actin polymerization force 
	calculateForce_AP(F_AP,X_t,x_t,c_B,mesh);	

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
	double area_tot = calculateAreaTot(x_t,mesh);
	Eigen::VectorXd F_areaConstraint(n_node*3); F_areaConstraint.setZero(); 
	calculateAreaConstraintForce(F_areaConstraint, X_t, x_t, mesh, area_tot);

	// net force
	Eigen::VectorXd f_t = F_bending + F_AP + F_FA + F_areaConstraint - Res;
	std::cout<<"net internal force norm "<<Res.norm()<<"\n";
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
			V_t[ni] = V_t_hlf[ni] + (dt_real - dt/2.0)*A_t[ni];
		}else{
			V_t[ni] = V_t[ni] + dt_real*A_t[ni];
		}
	}

}