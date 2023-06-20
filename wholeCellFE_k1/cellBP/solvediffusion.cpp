
#include "vec.h"
#include "nondestructivetrimesh.h"
#include <vector>
#include <Eigen/Dense> // most of the vector functions I will need inside of an element
#include <Eigen/Sparse> // functions for solution of linear systems
#include <Eigen/OrderingMethods>
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;
using namespace Eigen;

void advance_diffusion(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, std::vector<double> &c_B, std::vector<double> &c_C, std::vector<double> &c_BC, double dt)
{

//%% get parameters
int ne = mesh.num_triangles();// % get number of element
int nn = x.size(); //            % get total node number
// int nne = 3;//                   % get number of node per element 3 for triangle mesh // Not used

// --------------constant parameters-----------------
double DB = 30; //            % diffusion rate of ligand (BMP)    (microns^2*s^-1)*60s*m^-1
double DC = 20; //            % diffusion rate of ligand (Chd)    (microns^2*s^-1)*60s*m^-1
double DBC= 13; //            % diffusion rate of ligand (BC)    (microns^2*s^-1)*60s*m^-1

double j1 = 5/60; //          % production rate of BMP          nM*m^-1
double j2 = 5/60; //          % production rate of Chd          nM*m^-1

double k2   = 0.001; //      % binding rates for BMP ligand and Chordin          nM^-1*m^-1
double k_2  = 0.012; //      % unbinding rates for BMP ligand and Chordin        m^-1

double decB = 0.00025; //  % decay rate of Ligand (BMP)    nM*m^-1
double decC = 0.00025; //  % decay rate of Chd             nM*m^-1

double ks_bmp=1;
double ks_chd=1;

double tol=1e-3;          //  % tolerance for Newton method
double max_it=100;       // max iteration number


// -------------initial conditions-------------------
// initializing the residual vector, stiffness matrix, sparse linear solver
VectorXd RR(nn*3); 
VectorXd dU(nn*3);
SparseMatrix<double, ColMajor> KK(nn*3,nn*3);
std::vector<T> KK_triplets;
SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> >   solver;

//% iteration for  Newton-Raphson method
//l------------------oop for iteration-----------------------
 for (int it=0;it<max_it;it++)
{    
  // clear residual and stiffness for iteration
  KK_triplets.clear();
  RR.setZero();

  //gass value for the current concetration
  std::vector<double> c_B0 = c_B;
  std::vector<double> c_C0 = c_C;
  std::vector<double> c_BC0 = c_BC;

// --------------loop for element start
// this code is for solve the diffusion reaction equatio
  for(int i=0;i<ne;i++){
  // for triangle i   
// ---------------------mesh information---------
    int node1index = mesh.get_triangle(i)[0];
    int node2index = mesh.get_triangle(i)[1];
    int node3index = mesh.get_triangle(i)[2];

// get coordinates
// coords of node 1
    double node1x = x[node1index][0];
    double node1y = x[node1index][1];
    double node1z = x[node1index][2];
    // coords of node 2
    double node2x = x[node2index][0];
    double node2y = x[node2index][1];
    double node2z = x[node2index][2];
     // coords of node 3
    double node3x = x[node3index][0];
    double node3y = x[node3index][1];
    double node3z = x[node3index][2];
    // get concentrations
    // node 1 conc
    Vector3d Uc_BMP;
    Uc_BMP(0) = c_B0[node1index];
    Uc_BMP(1) = c_B0[node2index];
    Uc_BMP(2) = c_B0[node3index];
    
    Vector3d Uc_CHD;
    Uc_CHD(0) = c_C0[node1index];
    Uc_CHD(1) = c_C0[node2index];
    Uc_CHD(2) = c_C0[node3index];

    Vector3d Uc_BC;
    Uc_BC(0) = c_BC0[node1index];
    Uc_BC(1) = c_BC0[node2index];
    Uc_BC(2) = c_BC0[node3index];
    // node 2 conc

    Vector3d Uc_BMPt;
    Uc_BMPt(0) = c_B[node1index];
    Uc_BMPt(1) = c_B[node2index];
    Uc_BMPt(2) = c_B[node3index];
    
    Vector3d Uc_CHDt;
    Uc_CHDt(0) = c_C[node1index];
    Uc_CHDt(1) = c_C[node2index];
    Uc_CHDt(2) = c_C[node3index];

    Vector3d Uc_BCt;
    Uc_BCt(0) = c_BC[node1index];
    Uc_BCt(1) = c_BC[node2index];
    Uc_BCt(2) = c_BC[node3index];

//set source term for BMP
    Vector3d s_BMP;
    s_BMP(0)=ks_bmp*node1x*j1; //   %initialize source terms for bmp !!!!!! modify function by chang of x
    s_BMP(1)=ks_bmp*node2x*j1; //   %initialize source terms for bmp !!!!!! modify function by chang of x
    s_BMP(2)=ks_bmp*node3x*j1; //   %initialize source terms for bmp !!!!!! modify function by chang of x

//set source term for CHD
    Vector3d s_CHD;
    s_CHD(0)=ks_chd*node1x*j2; //   %initialize source terms for bmp !!!!!! modify function by chang of x
    s_CHD(1)=ks_chd*node2x*j2; //   %initialize source terms for bmp !!!!!! modify function by chang of x
    s_CHD(2)=ks_chd*node3x*j2; //   %initialize source terms for bmp !!!!!! modify function by chang of x
      
    Vector3d  dNxi;
    dNxi(0)=-1;
    dNxi(1)=1;
    dNxi(2)=0;
    Vector3d dNeta;
    dNeta(0)=-1;
    dNeta(1)=0;
    dNeta(2)=1;
    Vector3d N;
    N(0)=1/3.0;
    N(0)=1/3.0;
    N(1)=1/3.0;

    //double GradientN[2][3];
    MatrixXd GradientN(2,3);
    GradientN(0,0)=dNxi(0);
    GradientN(0,1)=dNxi(1);
    GradientN(0,2)=dNxi(2);
    GradientN(1,0)=dNeta(0);
    GradientN(1,1)=dNeta(1);
    GradientN(1,2)=dNeta(2);

    Vector3d g1,g2;
    for(int j=0;j<2;j++)
    {
      g1(j)=dNxi(j)*node1x+dNxi(j)*node1y+dNxi(j)*node1z;
      g2(j)=dNeta(j)*node1x+dNeta(j)*node1y+dNeta(j)*node1z;
    }


    //inverse matrix 
    Vector3d e1=g1/g1.norm();
    Vector3d g1xg2 = g1.cross(g2);
    Vector3d n = g1xg2/g1xg2.norm();
    Vector3d e2=n.cross(e1);
        
    Matrix3d A;
    A.col(0) = e1;
    A.col(1) = e2;
    A.col(2) = n;
    
    Matrix3d px;
    px(0,0)=node1x;
    px(1,0)=node1y;
    px(2,0)=node1z;
    px(0,1)=node2x;
    px(1,1)=node2y;
    px(2,1)=node2z;
    px(0,2)=node3x;
    px(1,2)=node3y;
    px(2,2)=node3z;

    // get new coordinate
    Matrix3d poiprime = A.inverse()*px;
    
    // took the first two row as the local coordinates for  surface diffusion
    MatrixXd poiprime2D(2,3);
    poiprime2D.row(0)=poiprime.row(0);
    poiprime2D.row(1)=poiprime.row(1);
    
    // get the element concentraion for current time point
    double c_BMP= N.dot(Uc_BMP);
    double c_CHD=N.dot(Uc_CHD);
    double c_BC=N.dot(Uc_BC);
    // get the element concentraion for test concentraion  
    double c_BMPt=N.dot(Uc_BMPt);
    double c_CHDt=N.dot(Uc_CHDt);
    double c_BCt=N.dot(Uc_BCt);
      
    Matrix2d J = GradientN*poiprime2D.transpose(); //% Jacobian Matrix
    double DetJ = J.determinant();//% Jacobian determinate
    MatrixXd B = J.inverse()*GradientN;
       
    Matrix2d eye = MatrixXd::Identity(2,2);
    Matrix2d D_BMP=DB*eye;//                 % apply diffusion rate
    Matrix2d D_CHD=DC*eye;//                 % apply diffusion rate
    Matrix2d D_BC=DBC*eye;//                 % apply diffusion rate

    Matrix3d ke_BMP=B.transpose()*D_BMP*B*DetJ;          //% element stiffness matrix
    Matrix3d ke_CHD=B.transpose()*D_CHD*B*DetJ;          //% element stiffness matrix
    Matrix3d ke_BC=B.transpose()*D_BC*B*DetJ;             //% element stiffness matrix
      
    double se_BMP=N.dot(s_BMP);//                % element source term of BMP
    double se_CHD=N.dot(s_CHD);//               % element source term of CHD
           
    //---------------------Set Residual-------------------------------------------------------------------------
    Vector3d Res_BME;
    Res_BME(0) = -N(0)*c_BMP/dt*DetJ+N(0)*c_BMPt/dt*DetJ-N(0)*(se_BMP-k2*c_BMPt*c_CHDt+k_2*c_BCt-decB*c_BMPt)*DetJ;
    Res_BME(1) = -N(1)*c_BMP/dt*DetJ+N(1)*c_BMPt/dt*DetJ-N(1)*(se_BMP-k2*c_BMPt*c_CHDt+k_2*c_BCt-decB*c_BMPt)*DetJ;
    Res_BME(2) = -N(2)*c_BMP/dt*DetJ+N(2)*c_BMPt/dt*DetJ-N(2)*(se_BMP-k2*c_BMPt*c_CHDt+k_2*c_BCt-decB*c_BMPt)*DetJ;
    Res_BME = Res_BME + ke_BMP*Uc_BMPt;

    Vector3d Res_CHD;
    Res_CHD(0) = -N(0)*c_CHD/dt*DetJ+N(0)*c_CHDt/dt*DetJ-N(0)*(se_CHD-k2*c_BMPt*c_CHDt+k_2*c_BCt-decC*c_CHDt)*DetJ;
    Res_CHD(1) = -N(1)*c_CHD/dt*DetJ+N(1)*c_CHDt/dt*DetJ-N(1)*(se_CHD-k2*c_BMPt*c_CHDt+k_2*c_BCt-decC*c_CHDt)*DetJ;
    Res_CHD(2) = -N(2)*c_CHD/dt*DetJ+N(2)*c_CHDt/dt*DetJ-N(2)*(se_CHD-k2*c_BMPt*c_CHDt+k_2*c_BCt-decC*c_CHDt)*DetJ;
    Res_CHD=Res_CHD+ke_CHD*Uc_CHDt;

    Vector3d Res_BC;
    Res_BC(0) = -N(0)*c_BC/dt*DetJ+N(0)*c_BCt/dt*DetJ-N(0)*(k2*c_BMPt*c_CHDt-k_2*c_BCt)*DetJ;
    Res_BC(1) = -N(1)*c_BC/dt*DetJ+N(1)*c_BCt/dt*DetJ-N(1)*(k2*c_BMPt*c_CHDt-k_2*c_BCt)*DetJ;
    Res_BC(2) = -N(2)*c_BC/dt*DetJ+N(2)*c_BCt/dt*DetJ-N(2)*(k2*c_BMPt*c_CHDt-k_2*c_BCt)*DetJ;
    Res_BC =Res_BC+ke_BC*Uc_BCt;

    Matrix3d K11;
    K11 << N(0)*N(0)*DetJ/dt-N(0)*(-k2*c_CHDt-decB)*N(0)*DetJ+ke_BMP(0,0), N(0)*N(1)*DetJ/dt-N(0)*(-k2*c_CHDt-decB)*N(1)*DetJ+ke_BMP(0,1), N(0)*N(2)*DetJ/dt-N(0)*(-k2*c_CHDt-decB)*N(2)*DetJ+ke_BMP(0,2),
           N(1)*N(0)*DetJ/dt-N(1)*(-k2*c_CHDt-decB)*N(0)*DetJ+ke_BMP(1,0), N(1)*N(1)*DetJ/dt-N(1)*(-k2*c_CHDt-decB)*N(1)*DetJ+ke_BMP(1,1), N(1)*N(2)*DetJ/dt-N(1)*(-k2*c_CHDt-decB)*N(2)*DetJ+ke_BMP(1,2),
           N(2)*N(0)*DetJ/dt-N(2)*(-k2*c_CHDt-decB)*N(0)*DetJ+ke_BMP(2,0), N(2)*N(1)*DetJ/dt-N(2)*(-k2*c_CHDt-decB)*N(1)*DetJ+ke_BMP(2,1), N(2)*N(2)*DetJ/dt-N(2)*(-k2*c_CHDt-decB)*N(2)*DetJ+ke_BMP(2,2);

    Matrix3d K22;  
    K22 << N(0)*N(0)*DetJ/dt-N(0)*(-k2*c_BMPt-decC)*N(0)*DetJ+ke_CHD(0,0), N(0)*N(1)*DetJ/dt-N(0)*(-k2*c_BMPt-decC)*N(1)*DetJ+ke_CHD(0,1), N(0)*N(2)*DetJ/dt-N(0)*(-k2*c_BMPt-decC)*N(2)*DetJ+ke_CHD(0,2),
           N(1)*N(0)*DetJ/dt-N(1)*(-k2*c_BMPt-decC)*N(0)*DetJ+ke_CHD(1,0), N(1)*N(1)*DetJ/dt-N(1)*(-k2*c_BMPt-decC)*N(1)*DetJ+ke_CHD(1,1), N(1)*N(2)*DetJ/dt-N(1)*(-k2*c_BMPt-decC)*N(2)*DetJ+ke_CHD(1,2),
           N(2)*N(0)*DetJ/dt-N(2)*(-k2*c_BMPt-decC)*N(0)*DetJ+ke_CHD(2,0), N(2)*N(1)*DetJ/dt-N(2)*(-k2*c_BMPt-decC)*N(1)*DetJ+ke_CHD(2,1), N(2)*N(2)*DetJ/dt-N(2)*(-k2*c_BMPt-decC)*N(2)*DetJ+ke_CHD(2,2);
    
    Matrix3d K33;   
    K33 << N(0)*N(0)*DetJ/dt-N(0)*(-k_2)*N(0)*DetJ+ke_BC(0,0), N(0)*N(1)*DetJ/dt-N(0)*(-k_2)*N(1)*DetJ+ke_BC(0,1), N(0)*N(2)*DetJ/dt-N(0)*(-k_2)*N(2)*DetJ+ke_BC(0,2),
           N(1)*N(0)*DetJ/dt-N(1)*(-k_2)*N(0)*DetJ+ke_BC(1,0), N(1)*N(1)*DetJ/dt-N(1)*(-k_2)*N(1)*DetJ+ke_BC(1,1), N(1)*N(2)*DetJ/dt-N(1)*(-k_2)*N(2)*DetJ+ke_BC(1,2),
           N(2)*N(0)*DetJ/dt-N(2)*(-k_2)*N(0)*DetJ+ke_BC(2,0), N(2)*N(1)*DetJ/dt-N(2)*(-k_2)*N(2)*DetJ+ke_BC(2,1), N(2)*N(2)*DetJ/dt-N(2)*(-k_2)*N(2)*DetJ+ke_BC(2,2);
    
    Matrix3d K12;  
    K12 << N(0)*(-k2*c_BMPt)*N(0)*DetJ,N(0)*(-k2*c_BMPt)*N(1)*DetJ,N(0)*(-k2*c_BMPt)*N(2)*DetJ,
           N(1)*(-k2*c_BMPt)*N(0)*DetJ,N(1)*(-k2*c_BMPt)*N(1)*DetJ,N(1)*(-k2*c_BMPt)*N(2)*DetJ,
           N(2)*(-k2*c_BMPt)*N(0)*DetJ,N(2)*(-k2*c_BMPt)*N(1)*DetJ,N(2)*(-k2*c_BMPt)*N(2)*DetJ;
    K12=-K12;

    Matrix3d K21;    
    K21 << N(0)*(-k2*c_CHDt)*N(0)*DetJ,N(0)*(-k2*c_CHDt)*N(1)*DetJ,N(0)*(-k2*c_CHDt)*N(2)*DetJ,
            N(1)*(-k2*c_CHDt)*N(0)*DetJ,N(1)*(-k2*c_CHDt)*N(1)*DetJ,N(1)*(-k2*c_CHDt)*N(2)*DetJ,
            N(2)*(-k2*c_CHDt)*N(0)*DetJ,N(2)*(-k2*c_CHDt)*N(1)*DetJ,N(2)*(-k2*c_CHDt)*N(2)*DetJ;
    K21=-K21;

    Matrix3d K13;    
    K13 << N(0)*(k_2)*N(0)*DetJ,N(0)*(k_2)*N(1)*DetJ,N(0)*(k_2)*N(2)*DetJ,
           N(1)*(k_2)*N(0)*DetJ,N(1)*(k_2)*N(1)*DetJ,N(1)*(k_2)*N(2)*DetJ,
           N(2)*(k_2)*N(0)*DetJ,N(2)*(k_2)*N(1)*DetJ,N(2)*(k_2)*N(2)*DetJ;
    K13=-K13;       
    
    Matrix3d K31;     
    K31 << N(0)*(k2*c_CHDt)*N(0)*DetJ,N(0)*(k2*c_CHDt)*N(1)*DetJ,N(0)*(k2*c_CHDt)*N(2)*DetJ,
           N(1)*(k2*c_CHDt)*N(0)*DetJ,N(1)*(k2*c_CHDt)*N(1)*DetJ,N(1)*(k2*c_CHDt)*N(2)*DetJ,
           N(2)*(k2*c_CHDt)*N(0)*DetJ,N(2)*(k2*c_CHDt)*N(1)*DetJ,N(2)*(k2*c_CHDt)*N(2)*DetJ; 
    K31=-K31; 
    
    Matrix3d K23; 
    K23=K13;

    Matrix3d K32;
    K32 << N(0)*(k2*c_BMPt)*N(0)*DetJ,N(0)*(k2*c_BMPt)*N(1)*DetJ,N(0)*(k2*c_BMPt)*N(2)*DetJ,
           N(1)*(k2*c_BMPt)*N(0)*DetJ,N(1)*(k2*c_BMPt)*N(1)*DetJ,N(1)*(k2*c_BMPt)*N(2)*DetJ,
           N(2)*(k2*c_BMPt)*N(0)*DetJ,N(2)*(k2*c_BMPt)*N(1)*DetJ,N(2)*(k2*c_BMPt)*N(2)*DetJ; 
    K32=-K32; 
      
    MatrixXd K_e(9,9); K_e.setZero();
    K_e.block<3,3>(0,0) = K11;
    K_e.block<3,3>(0,3) = K12;
    K_e.block<3,3>(0,6) = K13;
    K_e.block<3,3>(3,0) = K21;
    K_e.block<3,3>(3,3) = K22;
    K_e.block<3,3>(3,6) = K23;
    K_e.block<3,3>(6,0) = K31;
    K_e.block<3,3>(6,3) = K32;
    K_e.block<3,3>(6,6) = K33; 
    //std::cout<<K_e.determinant()<<"\n\n";

    VectorXd Res_e(9);Res_e.setZero();
    Res_e.head<3>() = Res_BME;
    Res_e.segment<3>(3) = Res_CHD;
    Res_e.tail<3>() = Res_BC;
      
// ASSEMBLY
    RR(node1index) += Res_e(0); // residual for cB for node 1 of element i
    RR(node2index) += Res_e(1); // residual for cB for node 2 of element i
    RR(node3index) += Res_e(2); // residual for cB for node 3 of element i

    RR(node1index+nn) += Res_e(3); // residual for cC for node 1 of element i
    RR(node2index+nn) += Res_e(4); // residual for cC for node 2 of element i
    RR(node3index+nn) += Res_e(5); // residual for cC for node 3 of element i

    RR(node1index+2*nn) += Res_e(6); // residual for cBC for node 1 of element i
    RR(node2index+2*nn) += Res_e(7); // residual for cBC for node 2 of element i
    RR(node3index+2*nn) += Res_e(8); // residual for cBC for node 3 of element i


    // KK_triplets.push_back(row,col,value)  //?????????????????????? necessary?
    //int varArray[] = {&node1index, &node2index, &node3index};
    std::vector<int> varArray(3,0);
    varArray[0] = node1index;
    varArray[1] = node2index;
    varArray[2] = node3index;

    for(int a=0;a<3;a++)
    {
      for(int b=0;b<3;b++)
      {
            KK_triplets.push_back(T(varArray[a],varArray[b],K_e(a,b)));
            KK_triplets.push_back(T(varArray[a],varArray[b]+nn,K_e(a,b+3)));
            KK_triplets.push_back(T(varArray[a],varArray[b]+2*nn,K_e(a,b+6)));
            KK_triplets.push_back(T(varArray[a]+nn,varArray[b],K_e(a+3,b)));
            KK_triplets.push_back(T(varArray[a]+nn,varArray[b]+nn,K_e(a+3,b+3)));
            KK_triplets.push_back(T(varArray[a]+nn,varArray[b]+2*nn,K_e(a+3,b+6)));
            KK_triplets.push_back(T(varArray[a]+2*nn,varArray[b],K_e(a+6,b)));
            KK_triplets.push_back(T(varArray[a]+2*nn,varArray[b]+nn,K_e(a+6,b+3)));
            KK_triplets.push_back(T(varArray[a]+2*nn,varArray[b]+2*nn,K_e(a+6,b+6)));
      }
          //KK_triplets.push_back(node1index,node1index,K_e(0,0));
          //KK_triplets.push_back(node1index,node2index,K_e(0,1));
          //KK_triplets.push_back(node1index,node3index,K_e(0,2));
          //KK_triplets.push_back(node2index,node1index,K_e(1,0));
          //KK_triplets.push_back(node2index,node2index,K_e(1,1));
          //KK_triplets.push_back(node2index,node3index,K_e(1,2));

          //KK_triplets.push_back(node1index,node1index+nn,K_e(0,3));
          //KK_triplets.push_back(node1index,node2index+nn,K_e(0,4));
          //KK_triplets.push_back(node1index,node3index+nn,K_e(0,5));

      //!!!!!!!!!!!!!!!!!!!!!!!!

    
     //need to know how to change variable name as in for loop


      //!!!!!!!!!!!!!!!!!!!!!!!!!
      //Kl(coni,coni)=K_e(1:3,1:3);           % gloable flag
      //Kl(coni+nn,coni+nn)=K_e(4:6,4:6);
      //Kl(coni+2*nn,coni+2*nn)=K_e(7:9,7:9);

      //K=K+Kl;
      //RES=RES+RESl;
    }
  }// closing the element loop
    double error;
    error=RR.norm();
    std::cout<<"Error = "<<error<<"\n";
    if(error>tol)
    {
      KK.setFromTriplets(KK_triplets.begin(), KK_triplets.end());
      KK.makeCompressed();
      //std::cout<<"KK\n"<<KK<<"\n\n";
      //std::cout<<"KK determinant\n"<<KK.determinant()<<"\n\n";
      //solver2.analyzePattern(KK2);
      // Compute the numerical factorization
      //solver2.factorize(KK2);
      solver.compute(KK);
      //Use the factors to solve the linear system
      dU = solver.solve(-1.*RR);
      //VectorXd dU(nn*3);
      //double dU=KK_triplets.inverse()*RR;
       //dU=-K\RES;
       //c_B0=c_B0+dU.head<nn>;
       //c_C0=c_B0+dU.segment<nn>(nn);
       //c_BC0=c_BC0+dU.tail<nn>;
      for(int i=0;i<nn;i++){
        c_B[i] += dU(i);
        c_C[i] += dU(i+nn);
        c_BC[i] += dU(i+2*nn);

      }
     }
    else
    {
       break;
    }
} // closes the Newton Raphson
} // closes the function 

