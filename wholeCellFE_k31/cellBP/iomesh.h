// ---------------------------------------------------------
//
//  iomesh.h
//
//  Non-member functions for reading and writing various mesh file formats.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_IOMESH_H
#define EL_TOPO_IOMESH_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <fstream>
#include <vec.h>
#include <vector>
#include <Eigen/Dense>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class NonDestructiveTriMesh;

namespace Gluvi
{
    struct Target3D;
}

// ---------------------------------------------------------
//  Function declarations
// ---------------------------------------------------------

// ---------------------------------------------------------
//
// Read/write mesh in our own binary format
//
// ---------------------------------------------------------

bool write_binary_file( const NonDestructiveTriMesh &mesh,  const std::vector<Vec3d> &x, const std::vector<double> &masses, double curr_t, const char *filename_format, ...);
bool write_binary_file_with_velocities( const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, const std::vector<double> &masses, const std::vector<Vec3d> &v, double curr_t, const char *filename_format, ...);

bool write_surface_ids( const std::vector<size_t> &ids,
                       const char *filename_format, ... );

bool read_binary_file( NonDestructiveTriMesh &mesh, std::vector<Vec3d> &x, std::vector<double> &masses, double& curr_t, const char *filename_format, ...);
bool read_binary_file_with_velocities( NonDestructiveTriMesh &mesh, std::vector<Vec3d> &x, std::vector<double> &masses, std::vector<Vec3d> &v, double& curr_t, const char *filename_format, ...);

bool read_surface_ids( std::vector<unsigned int> &ids,
                      const char *filename_format, ... );


// ---------------------------------------------------------
//
// Read/write mesh in Wavefront OBJ format (ASCII)
//
// ---------------------------------------------------------

bool write_objfile(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, const char *filename_format, ...);
bool read_objfile(NonDestructiveTriMesh &mesh, std::vector<Vec3d> &x, const char *filename_format, ...);

// ---------------------------------------------------------
//
// read (edited) COMSOL file
//
//----------------------------------------------------------
//

// by Adrian BT on July 2017
void read_ecomsol(NonDestructiveTriMesh &mesh, std::vector<Vec3d> &x, const char* filename);
void read_ecomsol(NonDestructiveTriMesh &mesh, std::vector<Vec3d> &x, std::vector<double> &c_B,std::vector<double> &c_C,std::vector<double> &c_BC, const char* filename);

// ---------------------------------------------------------
//
// write VTK file
//
//----------------------------------------------------------
//

// by Adrian BT on July 2017
void write_vtk(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, double current_time, const char* filename);
//void write_vtk(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, const std::vector<double> &c_B, const std::vector<double> &c_C, const std::vector<double> &c_BC, const std::vector<Vec3d> &fa_u, const std::vector<Eigen::VectorXd> &ea_t, double current_time, const char* filename, double kint);
void write_vtk(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, const std::vector<double> &c_A, const std::vector<double> &c_B, const std::vector<double> &c_C, const std::vector<Vec3d> &fa_u, const std::vector<Eigen::VectorXd> &ea_t, double current_time, const char* filename);
void write_vtk(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x,const std::vector<Eigen::VectorXd> &ea_t, const std::vector<Vec3d> &fa_u, double current_time, const char* filename);
// ---------------------------------------------------------
//
// Write mesh in Renderman RIB format (geometry only)
//
// ---------------------------------------------------------

bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, const char *filename_format, ...);
bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, std::ostream &output);
bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, FILE *output);

// ---------------------------------------------------------
///
/// Write an RIB file for the shadow map for the given light
///
// ---------------------------------------------------------

bool output_shadow_rib( Gluvi::Target3D& light, const std::vector<Vec3d>& positions,  const NonDestructiveTriMesh& mesh, const char *filename_format, ...);

// ---------------------------------------------------------
///
/// Write a render-ready RIB file.
///
// ---------------------------------------------------------

bool output_rib( const std::vector<Vec3d>& positions, const NonDestructiveTriMesh& mesh, const char *filename_format, ...);

// ---------------------------------------------------------
//
// Write mesh in PBRT format (geometry only)
//
// ---------------------------------------------------------

bool write_pbrtfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, const char *filename_format, ...);
bool write_pbrtfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, std::ostream &output);


// ---------------------------------------------------------
//
// Write an STL vector to an ASCII file.  Not really mesh-related, but useful.
//
// ---------------------------------------------------------

template<class T> void dump_vector_to_file( const char* filename, const std::vector<T, std::allocator<T> >& vec )
{
    std::ofstream outfile( filename, std::ios::out|std::ios::trunc );
    for ( unsigned int i = 0; i < vec.size(); ++i )
    {
        outfile << vec[i] << std::endl;
    }         
    outfile.close();
}


#endif
