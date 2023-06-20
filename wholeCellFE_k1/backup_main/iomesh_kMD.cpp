
// ---------------------------------------------------------
//
//  iomesh.cpp
//
//  Non-member functions for reading and writing various mesh file formats.
//
// ---------------------------------------------------------


#include <iomesh.h>

#include <bfstream.h>
#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <nondestructivetrimesh.h>



#ifndef NO_GUI
#include <gluvi.h>
#endif

#define LINESIZE 1024 // maximum line size when reading .OBJ files


// ---------------------------------------------------------
///
/// Write mesh in binary format
///
// ---------------------------------------------------------

bool write_binary_file( const NonDestructiveTriMesh &mesh, 
                       const std::vector<Vec3d> &x,
                       const std::vector<double> &masses, 
                       double curr_t, 
                       const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);   
    bofstream outfile( filename_format, ap );
    va_end(ap);
    
    outfile.write_endianity();
    
    outfile << curr_t;
    
    unsigned int nverts = static_cast<unsigned int>( x.size() );
    outfile << nverts;
    for ( unsigned int i = 0; i < x.size(); ++i )
    {
        outfile << x[i][0];
        outfile << x[i][1];
        outfile << x[i][2];
    }
    
    assert( x.size() == masses.size() );
    
    for ( unsigned int i = 0; i < masses.size(); ++i )
    {
        outfile << masses[i];
    }
    
    unsigned int ntris = static_cast<unsigned int>( mesh.num_triangles() );
    outfile << ntris;
    
    for ( unsigned int t = 0; t < mesh.num_triangles(); ++t )
    {
        const Vec3ui& tri = Vec3ui( mesh.get_triangle(t) );
        outfile << tri[0];
        outfile << tri[1];
        outfile << tri[2];      
    }
    
    outfile.close();
    
    return outfile.good();
}

// ---------------------------------------------------------
///
/// Write mesh in binary format, with per-vertex velocities
///
// ---------------------------------------------------------

bool write_binary_file_with_velocities( const NonDestructiveTriMesh &mesh, 
                                       const std::vector<Vec3d> &x,
                                       const std::vector<double> &masses,                                       
                                       const std::vector<Vec3d> &v,
                                       double curr_t, 
                                       const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);   
    bofstream outfile( filename_format, ap );
    va_end(ap);
    
    outfile.write_endianity();
    
    outfile << curr_t;
    
    unsigned int nverts = static_cast<unsigned int>( x.size() );
    outfile << nverts;
    
    for ( unsigned int i = 0; i < x.size(); ++i )
    {
        outfile << x[i][0];
        outfile << x[i][1];
        outfile << x[i][2];
    }
    
    assert( x.size() == masses.size() );
    
    for ( unsigned int i = 0; i < masses.size(); ++i )
    {
        outfile << masses[i];
    }
    
    for ( unsigned int i = 0; i < v.size(); ++i )
    {
        outfile << v[i][0];
        outfile << v[i][1];
        outfile << v[i][2];
    }
    
    unsigned int ntris = static_cast<unsigned int>( mesh.num_triangles() );
    outfile << ntris;
    
    for ( unsigned int t = 0; t < mesh.num_triangles(); ++t )
    {
        const Vec3ui& tri = Vec3ui( mesh.get_triangle(t) );
        outfile << tri[0];
        outfile << tri[1];
        outfile << tri[2];      
    }
    
    outfile.close();
    
    return outfile.good();
}

// ---------------------------------------------------------

bool write_surface_ids( const std::vector<size_t> &ids,
                       const char *filename_format, ... )
{
    va_list ap;
    va_start(ap, filename_format);   
    bofstream outfile( filename_format, ap );
    va_end(ap);
    
    outfile.write_endianity();
    
    size_t nids = ids.size(); 
    outfile << nids;
    
    for ( unsigned int i = 0; i < ids.size(); ++i )
    {
        unsigned int curr_id = static_cast<unsigned int>(ids[i]);
        outfile << curr_id;
    }
    
    outfile.close();
    
    return outfile.good();
}


// ---------------------------------------------------------
///
/// Read mesh in binary format
///
// ---------------------------------------------------------

bool read_binary_file( NonDestructiveTriMesh &mesh, 
                      std::vector<Vec3d> &x, 
                      std::vector<double> &masses, 
                      double& curr_t, 
                      const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);   
    bifstream infile( filename_format, ap );
    va_end(ap);
    
    assert( infile.good() );
    
    infile.read_endianity();
    
    infile >> curr_t;
    
    unsigned int nverts;
    infile >> nverts;
    x.resize( nverts );
    for ( unsigned int i = 0; i < nverts; ++i )
    {
        infile >> x[i][0];
        infile >> x[i][1];
        infile >> x[i][2];  
        mesh.add_vertex();
    }
    
    masses.resize( nverts );
    for ( unsigned int i = 0; i < nverts; ++i )
    {
        infile >> masses[i];
    }
    
    unsigned int ntris;
    infile >> ntris;
    
    for ( unsigned int t = 0; t < ntris; ++t )
    {
        Vec3ui tri;
        infile >> tri[0];
        infile >> tri[1];
        infile >> tri[2];
        
        mesh.add_triangle( Vec3st(tri) );
    }
    
    infile.close();
    
    return infile.good();
}


// ---------------------------------------------------------
///
/// Read mesh in binary format, with per-vertex velocities
///
// ---------------------------------------------------------

bool read_binary_file_with_velocities( NonDestructiveTriMesh &mesh, 
                                      std::vector<Vec3d> &x, 
                                      std::vector<double> &masses,
                                      std::vector<Vec3d> &v, 
                                      double& curr_t, 
                                      const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);   
    bifstream infile( filename_format, ap );
    va_end(ap);
    
    infile.read_endianity();
    
    infile >> curr_t;
    
    unsigned int nverts;
    infile >> nverts;
    
    mesh.set_num_vertices(nverts);
    
    x.resize( nverts );
    for ( unsigned int i = 0; i < nverts; ++i )
    {
        infile >> x[i][0];
        infile >> x[i][1];
        infile >> x[i][2];      
    }
    
    masses.resize( nverts );
    for ( unsigned int i = 0; i < masses.size(); ++i )
    {
        infile >> masses[i];
    }
    
    v.resize( nverts );
    for ( unsigned int i = 0; i < nverts; ++i )
    {
        infile >> v[i][0];
        infile >> v[i][1];
        infile >> v[i][2];      
    }
    
    
    unsigned int ntris;
    infile >> ntris;
    
    for ( unsigned int t = 0; t < ntris; ++t )
    {
        Vec3st tri;
        infile >> tri[0];
        infile >> tri[1];
        infile >> tri[2];
        mesh.add_triangle(tri);
    }
    
    infile.close();
    
    return infile.good();
}


// ---------------------------------------------------------

bool read_surface_ids( std::vector<unsigned int>& ids,
                      const char *filename_format, ... )
{
    
    va_list ap;
    va_start(ap, filename_format);   
    bifstream infile( filename_format, ap );
    va_end(ap);
    
    infile.read_endianity();
    
    unsigned int n;
    infile >> n;
    
    ids.resize(n);
    for ( unsigned int i = 0; i < ids.size(); ++i )
    {
        infile >> ids[i];
    }
    
    infile.close();
    
    return infile.good();
    
}


// ---------------------------------------------------------
///
/// Write mesh in Wavefront OBJ format
///
// ---------------------------------------------------------

bool write_objfile(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, const char *filename_format, ...)
{

#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
    +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);
    std::ofstream output(filename, std::ofstream::binary);

    delete[] filename;
    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);

    std::ofstream output(filename, std::ofstream::binary);
    std::free(filename);
    va_end(ap);
#endif    
    if(!output.good()) return false;
    
    output<<"# generated by editmesh"<<std::endl;
    for(unsigned int i=0; i<x.size(); ++i)
        output<<"v "<<x[i]<<std::endl;
    for(unsigned int t=0; t<mesh.num_triangles(); ++t)
        output<<"f "<<mesh.get_triangle(t)[0]+1<<' '<<mesh.get_triangle(t)[1]+1<<' '<<mesh.get_triangle(t)[2]+1<<std::endl; // correct for 1-based indexing in OBJ files
    return output.good();
}

namespace {
    
    // ---------------------------------------------------------
    ///
    /// Helper for reading OBJ file
    ///
    // ---------------------------------------------------------
    
    bool read_int(const char *s, int &value, bool &leading_slash, int &position)
    {
        leading_slash=false;
        for(position=0; s[position]!=0; ++position){
            switch(s[position]){
                case '/':
                    leading_slash=true;
                    break;
                case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
                    goto found_int;
            }
        }
        return false;
        
    found_int:
        value=0;
        for(;; ++position){
            switch(s[position]){
                case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
                    value=10*value+s[position]-'0';
                    break;
                default:
                    return true;
            }
        }
        return true; // should never get here, but keeps compiler happy
    }
    
    // ---------------------------------------------------------
    ///
    /// Helper for reading OBJ file
    ///
    // ---------------------------------------------------------
    
    void read_face_list(const char *s, std::vector<int> &vertex_list)
    {
        vertex_list.clear();
        int v, skip;
        bool leading_slash;
        for(int i=0;;){
            if(read_int(s+i, v, leading_slash, skip)){
                if(!leading_slash)
                    vertex_list.push_back(v-1); // correct for 1-based index
                i+=skip;
            }else
                break;
        }
    }
    
}  // unnamed namespace

// ---------------------------------------------------------
///
/// Read mesh in Wavefront OBJ format
///
// ---------------------------------------------------------

bool read_objfile(NonDestructiveTriMesh &mesh, std::vector<Vec3d> &x, const char *filename_format, ...)
{

#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);
    std::ifstream input(filename, std::ifstream::binary);
    delete[] filename;
    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);
    std::ifstream input(filename, std::ifstream::binary);
    std::free(filename);
    va_end(ap);
#endif

    if(!input.good()) return false;
    
    x.clear();
    mesh.clear();
    
    char line[LINESIZE];
    std::vector<int> vertex_list;
    while(input.good()){
        input.getline(line, LINESIZE);
        switch(line[0]){
            case 'v': // vertex data
                if(line[1]==' '){
                    Vec3d new_vertex;
                    std::sscanf(line+2, "%lf %lf %lf", &new_vertex[0], &new_vertex[1], &new_vertex[2]);
                    x.push_back(new_vertex);
                }
                break;
            case 'f': // face data
                if(line[1]==' '){
                    read_face_list(line+2, vertex_list);
                    for(int j=0; j<(int)vertex_list.size()-2; ++j)
                        mesh.add_triangle( Vec3st(vertex_list[0], vertex_list[j+1], vertex_list[j+2]) );
                }
                break;
        }
    }
    return true;
}

// ---------------------------------------------------------
///
/// Read mesh from (Edited) COMSOL file
///
// ---------------------------------------------------------

void read_ecomsol(NonDestructiveTriMesh &mesh, std::vector<Vec3d> &x, const char* filename)
{
    // read nodes
    x.clear();
    std::ifstream myfile(filename);
    std::string line;
    std::string keyword_node = "Nodes";
    int n_node;
    if (myfile.is_open())
    {
        // read in until you find the keyword Nodes
        while ( getline (myfile,line) )
        {
            // check for the keyword
            std::size_t found = line.find(keyword_node);
            if (found!=std::string::npos)
            {
                getline (myfile,line);
                std::stringstream ss0(line);
                ss0>>n_node;
                // found the beginning of the nodes, read the nnumber of nodes
                for(int i =0;i<n_node;i++){
                    getline (myfile,line);
                    std::stringstream ss1(line); 
                    Vec3d new_vertex;
                    ss1>>new_vertex[0];
                    ss1>>new_vertex[1];
                    ss1>>new_vertex[2];
                    new_vertex[2]=0; // for planar problems 
                    x.push_back(new_vertex);
                }
            }
        }
    }
    myfile.close();

    // read triangles
    int n_tri;
    myfile.open(filename);
    std::string keyword_element = "Elements";
    int auxE0,auxE1,auxE2;
    if (myfile.is_open())
    {
        // read in until you find the keyword Elements
        while ( getline (myfile,line) )
        {
            // check for the keyword
            std::size_t found = line.find(keyword_element);
            if (found!=std::string::npos)
            {
                // found the beginning of the elements, get number of elements and then loop
                getline (myfile,line);
                std::stringstream ss3(line);
                ss3>>n_tri;
                for(int i=0;i<n_tri;i++)
                {
                    getline (myfile,line);
                    std::stringstream ss4(line);
                    ss4>>auxE0;
                    ss4>>auxE1;
                    ss4>>auxE2;
                    mesh.add_triangle( Vec3st(auxE0, auxE1, auxE2) );
                }
            }
        }
    }
    myfile.close();

    // just to check, write it immediately as a VTK
    std::string filename2(filename);
    filename2.erase(filename2.size()-3);
    std::string filename3 = filename2+"vtk";
    write_vtk(mesh,x,-1,filename3.c_str());
    
}

// ---------------------------------------------------------
///
/// Read mesh from (Edited) COMSOL file plus concentrations
///
// ---------------------------------------------------------

void read_ecomsol(NonDestructiveTriMesh &mesh, std::vector<Vec3d> &x, std::vector<double> &c_B,std::vector<double> &c_C,std::vector<double> &c_BC, const char* filename)
{
    // read nodes
    x.clear();
    std::ifstream myfile(filename);
    std::string line;
    std::string keyword_node = "Nodes";
    int n_node;
    if (myfile.is_open())
    {
        // read in until you find the keyword Nodes
        while ( getline (myfile,line) )
        {
            // check for the keyword
            std::size_t found = line.find(keyword_node);
            if (found!=std::string::npos)
            {
                getline (myfile,line);
                std::stringstream ss0(line);
                ss0>>n_node;
                // found the beginning of the nodes, read the nnumber of nodes
                for(int i =0;i<n_node;i++){
                    getline (myfile,line);
                    std::stringstream ss1(line); 
                    Vec3d new_vertex;
                    ss1>>new_vertex[0];
                    ss1>>new_vertex[1];
                    ss1>>new_vertex[2];
                    x.push_back(new_vertex);
                }
            }
        }
    }
    myfile.close();

    // read triangles
    int n_tri;
    myfile.open(filename);
    std::string keyword_element = "Elements";
    int auxE0,auxE1,auxE2;
    if (myfile.is_open())
    {
        // read in until you find the keyword Elements
        while ( getline (myfile,line) )
        {
            // check for the keyword
            std::size_t found = line.find(keyword_element);
            if (found!=std::string::npos)
            {
                // found the beginning of the elements, get number of elements and then loop
                getline (myfile,line);
                std::stringstream ss3(line);
                ss3>>n_tri;
                for(int i=0;i<n_tri;i++)
                {
                    getline (myfile,line);
                    std::stringstream ss4(line);
                    ss4>>auxE0;
                    ss4>>auxE1;
                    ss4>>auxE2;
                    mesh.add_triangle( Vec3st(auxE0, auxE1, auxE2) );
                }
            }
        }
    }
    myfile.close();

    // read concentrations
    c_B.clear();
    c_C.clear();
    c_BC.clear();
    myfile.open(filename);
    std::string keyword_concentrations = "Concentrations";
    int n_conc;
    if (myfile.is_open())
    {
        // read in until you find the keyword 
        while ( getline (myfile,line) )
        {
            // check for the keyword
            std::size_t found = line.find(keyword_concentrations);
            if (found!=std::string::npos)
            {
                getline (myfile,line);
                std::stringstream ss0(line);
                ss0>>n_conc;
                // found the beginning of the nodal concentrentrations, read
                for(int i =0;i<n_conc;i++){
                    getline (myfile,line);
                    std::stringstream ss1(line); 
                    double cBi,cCi,cBCi;
                    ss1>>cBi;
                    ss1>>cCi;
                    ss1>>cBCi;
                    c_B.push_back(cBi);
                    c_C.push_back(cCi);
                    c_BC.push_back(cBCi);
                }
            }
        }
    }
    myfile.close();

    // just to check, write it immediately as a VTK
    std::string filename2(filename);
    filename2.erase(filename2.size()-3);
    std::string filename3 = filename2+"vtk";
    write_vtk(mesh,x,-1,filename3.c_str());
    
}

// ---------------------------------------------------------
///
/// Write VTK file
///
// ---------------------------------------------------------

void write_vtk(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, double current_time, const char* filename)
{
    std::ofstream savefile(filename);
    if (!savefile) {
        throw std::runtime_error("Unable to open output file.");
    }
    savefile<<"# vtk DataFile Version 2.0\nCellBP Time: "<<current_time<<"\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    savefile<<"POINTS "<<x.size()<<" double\n";
    for(int i=0;i<x.size();i++)
    {
        savefile<<x[i][0]<<" "<<x[i][1]<<" "<<x[i][2]<<"\n";
    }
    
    savefile<<"CELLS "<<mesh.num_triangles()<<" "<<mesh.num_triangles()*4<<"\n";
    for(int i=0;i<mesh.num_triangles();i++)
    {
        savefile<<"3";
        for(int j=0;j<3;j++)
        {
            savefile<<" "<<mesh.get_triangle(i)[j];
        }
        savefile<<"\n";
    }
    savefile<<"CELL_TYPES "<<mesh.num_triangles()<<"\n";
    for(int i=0;i<mesh.num_triangles();i++)
    {
        savefile<<"5\n";
    }


    
    savefile.close();
}

// ---------------------------------------------------------
///
/// Write VTK file with strains
///
// ---------------------------------------------------------

void write_vtk(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x,const std::vector<Eigen::VectorXd> &ea_t, const std::vector<Vec3d> &fa_u, double current_time, const char* filename)
{
    std::ofstream savefile(filename);
    if (!savefile) {
        throw std::runtime_error("Unable to open output file.");
    }
    savefile<<"# vtk DataFile Version 2.0\nCellBP Time: "<<current_time<<"\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    savefile<<"POINTS "<<x.size()<<" double\n";
    savefile<<std::fixed;
    for(int i=0;i<x.size();i++)
    {
        savefile<<x[i][0]<<" "<<x[i][1]<<" "<<x[i][2]<<"\n";
    }
    
    savefile<<"CELLS "<<mesh.num_triangles()<<" "<<mesh.num_triangles()*4<<"\n";
    for(int i=0;i<mesh.num_triangles();i++)
    {
        savefile<<"3";
        for(int j=0;j<3;j++)
        {
            savefile<<" "<<mesh.get_triangle(i)[j];
        }
        savefile<<"\n";
    }
    savefile<<"CELL_TYPES "<<mesh.num_triangles()<<"\n";
    for(int i=0;i<mesh.num_triangles();i++)
    {
        savefile<<"5\n";
    }

    savefile<<"POINT_DATA "<<x.size()<<"\n";
    savefile<<"SCALARS theta_e float 1\n";
    savefile<<"LOOKUP_TABLE default\n";
    for(int i=0;i<x.size();i++)
    {
        // Get normal
        Eigen::Vector3d normal; normal.setZero();
        // Get the triangles incident in this vertex
        double Asum = 0.0;
        std::vector<size_t> onering = mesh.m_vertex_to_triangle_map[i];
        for(int ti=0;ti<onering.size();ti++){
            // get the normal and area of each tringle
            const Vec3st& tri = mesh.get_triangle(onering[ti]);
            Eigen::Vector3d tri_n0; tri_n0<<x[tri[0]][0],x[tri[0]][1],x[tri[0]][2] ;
            Eigen::Vector3d tri_n1; tri_n1<<x[tri[1]][0],x[tri[1]][1],x[tri[1]][2] ;
            Eigen::Vector3d tri_n2; tri_n2<<x[tri[2]][0],x[tri[2]][1],x[tri[2]][2] ;
            Eigen::Vector3d E1 = tri_n1 - tri_n0; 
            Eigen::Vector3d E2 = tri_n2 - tri_n0; 
            Eigen::Vector3d normal_ti = E1.cross(E2);
            double area_ti = 0.5*normal_ti.norm();
            normal_ti = normal_ti/normal_ti.norm();
            
            // add to the average calculation
            normal += area_ti*normal_ti;
            Asum += area_ti/3; 
        }
        normal = normal/Asum;
        Eigen::Vector3d N= normal/normal.norm();
        
        Eigen::Matrix3d ea_t_m; 
        ea_t_m<<ea_t[i](0),ea_t[i](3),ea_t[i](4),
                ea_t[i](3),ea_t[i](1),ea_t[i](5),
                ea_t[i](4),ea_t[i](5),ea_t[i](2);

        // will need to go from ea_t to be_t to get the area change
        Eigen::Matrix3d Id; Id<<1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0;
        Eigen::Matrix3d aux = Id-2.0*ea_t_m;
        Eigen::Matrix3d be_t = aux.inverse();
        Eigen::Matrix3d NdyadN = N*N.transpose();
        Eigen::Matrix3d be_t_s = be_t - NdyadN ;
    
        // Area changes to get the current normal strain 
        double theta_e_t = sqrt(be_t.determinant()); // use the total be_t with 1*NdyadN
    
        savefile<<theta_e_t<<"\n";
    }
    savefile<<"VECTORS fa_u float\n";
    // for(int ni=0;ni<x.size();ni++){
    //     Eigen::Vector3d normal; normal.setZero();
    //     // Get the triangles inciden in this vertex
    //     double Asum = 0.0;
    //     std::vector<size_t> onering = mesh.m_vertex_to_triangle_map[ni];
    //     for(int ti=0;ti<onering.size();ti++){
    //         // get the normal and area of each tringle
    //         const Vec3st& tri = mesh.get_triangle(onering[ti]);
    //         Eigen::Vector3d tri_n0; tri_n0<<x[tri[0]][0],x[tri[0]][1],x[tri[0]][2] ;
    //         Eigen::Vector3d tri_n1; tri_n1<<x[tri[1]][0],x[tri[1]][1],x[tri[1]][2] ;
    //         Eigen::Vector3d tri_n2; tri_n2<<x[tri[2]][0],x[tri[2]][1],x[tri[2]][2] ;
    //         Eigen::Vector3d E1 = tri_n1 - tri_n0; 
    //         Eigen::Vector3d E2 = tri_n2 - tri_n0; 
    //         Eigen::Vector3d normal_ti = E1.cross(E2);
    //         double area_ti = 0.5*normal_ti.norm();
    //         normal_ti = normal_ti/normal_ti.norm();
            
    //         // add to the average calculation
    //         normal += area_ti*normal_ti;
    //         Asum += area_ti; 
    //     }
    //     normal = normal/Asum;
    //     savefile<<normal(0)<<" "<<normal(1)<<" "<<normal(2)<<"\n";
    // }
    for(int ni=0;ni<x.size();ni++){
        savefile<<fa_u[ni][0]<<" "<<fa_u[ni][1]<<" "<<fa_u[ni][2]<<"\n";
    }
    
    savefile.close();
}

// ---------------------------------------------------------
///
/// Write VTK file with concentrations
///
// ---------------------------------------------------------

void write_vtk(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, const std::vector<double> &c_A,const std::vector<double> &c_B,const std::vector<double> &c_C,const std::vector<Vec3d> &fa_u, const std::vector<Eigen::VectorXd> &ea_t, double current_time, const char* filename, double kint)
{
    std::ofstream savefile(filename);
    if (!savefile) {
        throw std::runtime_error("Unable to open output file.");
    }
    savefile<<"# vtk DataFile Version 2.0\nCellBP Time: "<<current_time<<"\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    savefile<<"POINTS "<<x.size()<<" double\n";
    for(int i=0;i<x.size();i++)
    {
        savefile<<x[i][0]<<" "<<x[i][1]<<" "<<x[i][2]<<"\n";
    }
    
    savefile<<"CELLS "<<mesh.num_triangles()<<" "<<mesh.num_triangles()*4<<"\n";
    for(int i=0;i<mesh.num_triangles();i++)
    {
        savefile<<"3";
        for(int j=0;j<3;j++)
        {
            savefile<<" "<<mesh.get_triangle(i)[j];
        }
        savefile<<"\n";
    }
    savefile<<"CELL_TYPES "<<mesh.num_triangles()<<"\n";
    for(int i=0;i<mesh.num_triangles();i++)
    {
        savefile<<"5\n";
    }

    savefile<<"POINT_DATA "<<c_B.size()<<"\n";
	savefile<<"SCALARS concentrations float 4\n";
	savefile<<"LOOKUP_TABLE default\n";
    //double k_FA = 1000; // pN/um
    double k_FA = 31000;
    double rho_i_max = 100; // max integrin density in #/um^2 from "Influence of type I collagen surface density on fibroblast spreading, motility, and contractility.
	for(int i=0;i<c_A.size();i++)
    {
        // Get normal
        Eigen::Vector3d normal; normal.setZero();
        // Get the triangles incident in this vertex
        double Asum = 0.0;
        std::vector<size_t> onering = mesh.m_vertex_to_triangle_map[i];
        for(int ti=0;ti<onering.size();ti++){
            // get the normal and area of each tringle
            const Vec3st& tri = mesh.get_triangle(onering[ti]);
            Eigen::Vector3d tri_n0; tri_n0<<x[tri[0]][0],x[tri[0]][1],x[tri[0]][2] ;
            Eigen::Vector3d tri_n1; tri_n1<<x[tri[1]][0],x[tri[1]][1],x[tri[1]][2] ;
            Eigen::Vector3d tri_n2; tri_n2<<x[tri[2]][0],x[tri[2]][1],x[tri[2]][2] ;
            Eigen::Vector3d E1 = tri_n1 - tri_n0; 
            Eigen::Vector3d E2 = tri_n2 - tri_n0; 
            Eigen::Vector3d normal_ti = E1.cross(E2);
            double area_ti = 0.5*normal_ti.norm();
            normal_ti = normal_ti/normal_ti.norm();
            
            // add to the average calculation
            normal += area_ti*normal_ti;
            Asum += area_ti/3; 
        }
        normal = normal/Asum;
        Eigen::Vector3d N= normal/normal.norm();
        
        Eigen::Matrix3d ea_t_m; 
        ea_t_m<<ea_t[i](0),ea_t[i](3),ea_t[i](4),
                ea_t[i](3),ea_t[i](1),ea_t[i](5),
                ea_t[i](4),ea_t[i](5),ea_t[i](2);

        // will need to go from ea_t to be_t to get the area change
        Eigen::Matrix3d Id; Id<<1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0;
        Eigen::Matrix3d aux = Id-2.0*ea_t_m;
        Eigen::Matrix3d be_t = aux.inverse();
        Eigen::Matrix3d NdyadN = N*N.transpose();
        Eigen::Matrix3d be_t_s = be_t - NdyadN ;
    
        // Area changes to get the current normal strain 
        // Careful with diffusion or not
        double theta_e_t = sqrt(be_t.determinant()); // use the total be_t with 1*NdyadN
        double FA_force_x = c_C[i]*k_FA*fa_u[i][0]*Asum*rho_i_max;
        double FA_force_y = c_C[i]*k_FA*fa_u[i][1]*Asum*rho_i_max;
        double normFA = sqrt(FA_force_x*FA_force_x + FA_force_y*FA_force_y );

        double force_s; // spring force vector
		double D  = sqrt(fa_u[i][0]*fa_u[i][0] + fa_u[i][1]*fa_u[i][1] + fa_u[i][2]*fa_u[i][2]); // displacement of the spring
		
		int fFA_npts = 6; // for 1 nm/ns pull rate piecewise linear fit
			
		Eigen::VectorXd pxFA(fFA_npts); pxFA.setZero();
		pxFA<< 0.,  2.51250583e-03,  4.40265201e-03,  5.91408228e-03, 11.76817302e-03, 17.001e-03; // extension in um
		
		Eigen::VectorXd pyFA(fFA_npts); pyFA.setZero();
		pyFA<< -7.15514497,  71.80501308, 356.62420075, 356.96656556, 126.21811585, 528.40481124; // force in pN

		Eigen::VectorXd k_MD(fFA_npts); k_MD.setZero();
		k_MD<< 31426.8556 , 99485.86235, 67080.67086, -7908.39043, 21978.66923, 76858.39733; // integrin stiffness in pN/um
		
		if(kint<0){
			for(int j=0;j<fFA_npts-1;j++){ // depending on the integrin displacement u_x, the spring stiffness is assigned
				if(D<pxFA[j+1] && D>=pxFA[j]){
					force_s = pyFA[j]+k_MD[j]*(D-pxFA[j]);
				}else if(D>=pxFA[fFA_npts-1]){
					force_s = pyFA[fFA_npts-1]+k_MD[fFA_npts-1]*(D-pxFA[fFA_npts-1]);
				}
			}
		}
		else{
			force_s = kint*D;
		}

        // double Fint_force_x = k_FA*fa_u[i][0];
        // double Fint_force_y = k_FA*fa_u[i][1];
        // double normFint = sqrt(Fint_force_x*Fint_force_x + Fint_force_y*Fint_force_y );

        // savefile<<c_A[i]<<" "<<theta_e_t<<" "<<c_C[i]<<" "<<normFint<<"\n";
        savefile<<c_A[i]<<" "<<theta_e_t<<" "<<c_C[i]<<" "<<force_s<<"\n";
    }
    savefile<<"VECTORS fa_u float\n";
    for(int ni=0;ni<x.size();ni++){
        // debug the normal, so just the boundary 
        savefile<<fa_u[ni][0]<<" "<<fa_u[ni][1]<<" "<<fa_u[ni][2]<<"\n";
        /*
        // if the boundary then let's save as the normal 
        // check if this vertex is on the boundary
        if(mesh.m_is_boundary_vertex[ni]){
            // get the vector in Eigen format 
            Eigen::Vector3d ni_x; ni_x<<x[ni][0],x[ni][1],x[ni][2];
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
                        Eigen::Vector3d pb; pb<<x[nj][0],x[nj][1],x[nj][2];
                        // actually need to know which is the interior node
                        int nint;
                        for(int enk=0;enk<3;enk++){
                            if(tri[enk]!=nj && tri[enk]!=ni){nint=tri[enk];break;}
                        }
                        Eigen::Vector3d pint; pint<<x[nint][0],x[nint][1],x[nint][2];
                        Eigen::Vector3d ncheck = (ni_x-pb).cross(pint-ni_x);
                        if(ncheck[2]>0){
                            ni_j1_x << x[nj][0],x[nj][1],x[nj][2];
                            n1check = 1;
                        }else{
                            ni_j2_x << x[nj][0],x[nj][1],x[nj][2];
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

            //savefile<<kappa*normal_i[0]<<" "<<kappa*normal_i[1]<<" "<<kappa*normal_i[2]<<"\n";
            
        }else{
            //savefile<<0.00000001<<" "<<0<<" "<<0<<"\n";
        }
        */
    }
    savefile.close();
}

// ---------------------------------------------------------
///
/// Write mesh in Renderman RIB format.
///
// ---------------------------------------------------------

bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, const char *filename_format, ...)
{
#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);
    std::ofstream output(filename, std::ofstream::binary);
    delete[] filename;
    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);
    std::ofstream output(filename, std::ofstream::binary);
    std::free(filename);
    va_end(ap);
#endif
    if(!output.good()) return false;
    return write_ribfile(mesh, x, output);
}

// ---------------------------------------------------------
///
/// Write mesh in Renderman RIB format.
///
// ---------------------------------------------------------

bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, std::ostream &output)
{
    output<<"# generated by editmesh"<<std::endl;
    output<<"PointsPolygons"<<std::endl;
    output<<" [ ";
    const std::vector<Vec3st>& tris = mesh.get_triangles();
    for(unsigned int i=0; i<tris.size(); ++i){
        output<<"3 ";
        if(i%38==37 && i!=tris.size()-1) output<<std::endl;
    }
    output<<"]"<<std::endl;
    output<<" [ ";
    for(unsigned int i=0; i<tris.size(); ++i){
        output<<tris[i]<<"  ";
        if(i%6==5 && i!=tris.size()-1) output<<std::endl;
    }
    output<<"]"<<std::endl;
    output<<" \"P\" [";
    for(unsigned int i=0; i<x.size(); ++i){
        output<<x[i]<<"  ";
        if(i%4==3 && i!=x.size()-1) output<<std::endl;
    }
    output<<"]"<<std::endl;
    
    return output.good();
}

// ---------------------------------------------------------
///
/// Write mesh in Renderman RIB format.
///
// ---------------------------------------------------------

bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, FILE *output)
{
    fprintf( output, "# generated by editmesh\n" );
    fprintf( output, "PointsPolygons\n" );
    fprintf( output, " [ " );
    const std::vector<Vec3st>& tris = mesh.get_triangles();
    for(unsigned int i=0; i<tris.size(); ++i){
        fprintf( output, "3 " );
        if(i%38==37 && i!=tris.size()-1) fprintf( output, "\n" );
    }
    fprintf( output, "]\n" );
    fprintf( output, " [ " );
    for(unsigned int i=0; i<tris.size(); ++i){
        Vec3ui new_tri;
        fprintf( output, " %d %d %d ", (int)tris[i][0], (int)tris[i][1], (int)tris[i][2] );
        if(i%6==5 && i!=tris.size()-1) fprintf( output, "\n" ); 
    }
    fprintf( output, "]\n" );
    fprintf( output, " \"P\" [" );
    for(unsigned int i=0; i<x.size(); ++i){
        fprintf( output, " %f ", x[i] );
        if(i%4==3 && i!=x.size()-1) fprintf( output, "\n" ); 
    }
    fprintf( output, "]\n" );
    
    return true; 
}


#ifndef NO_GUI

// ---------------------------------------------------------
///
/// Write an RIB file for the shadow map for the given light
///
// ---------------------------------------------------------

bool output_shadow_rib( Gluvi::Target3D& light, const std::vector<Vec3d>& positions, const NonDestructiveTriMesh& mesh, const char *filename_format, ...)
{
#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);

    std::ofstream out;
    out.open( filename );

    delete[] filename;

    if( out == NULL )
    {
        return false;
    }

    len=_vscprintf("track%04d_shadow.tiff", ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    filename=new char[len];
    vsprintf( filename, "track%04d_shadow.tiff", ap );

    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);   
    
    std::ofstream out;
    out.open( filename );
    
    delete[] filename;
    
    if( out == NULL )
    {
        return false;
    }
    
    vasprintf( &filename, "track%04d_shadow.tiff", ap );
    
    va_end(ap);
#endif
    // flatten
    std::vector<float> xs;
    for ( unsigned int i = 0; i < positions.size(); ++i )
    {
        xs.push_back( (float) positions[i][0] );
        xs.push_back( (float) positions[i][1] );
        xs.push_back( (float) positions[i][2] );
    }
    
    out << "Display \"" << filename << "\" \"file\" \"z\"" << std::endl;
    delete[] filename;
    
    // next line: image format (width and height in pixels, pixel aspect ratio)
    out << "Format " << "1024 1024 1" << std::endl;
    out << "PixelFilter \"box\" 1 1 " << std::endl;
    
    // then write out the camera specification
    light.export_rib(out);
    
    // start the scene
    out << "WorldBegin\n";
    
    out << "AttributeBegin\n";
    out << "  Color [0.6 0.6 0.2]\n";
    out << "  Opacity [1 1 1]\n";
    out << "  Surface \"matte\" \"Kd\" 1\n";
    
    const float plane_limit = 50.0f;
    char buf[256];
    sprintf( buf, "  Polygon \"P\" [-%f -%f -%f %f -%f -%f  %f %f -%f  -%f %f -%f]\n", 
            plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit );
    out << buf;
    
    out << "AttributeEnd\n";
    
    out << "Color [0.7 0.7 0.9]\n";
    out << "Surface \"matte\" \n";
    
    write_ribfile( mesh, xs, out);
    
    // finish the scene
    out << "WorldEnd\n";
    
    out.flush();
    out.close();
    
    return true;
}   


// ---------------------------------------------------------
///
/// Write a render-ready RIB file.
///
// ---------------------------------------------------------

bool output_rib( const std::vector<Vec3d>& positions, const NonDestructiveTriMesh& mesh, const char *filename_format, ...)
{
#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);

    std::ofstream out;
    out.open( filename );

    delete[] filename;

    if( out == NULL )
    {
        return false;
    }

    // first line: what image file this RIB file should produce
    len=_vscprintf("track%04d.tiff", ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    filename=new char[len];
    vsprintf( filename, "track%04d.tiff", ap );

    len=_vscprintf("track%04d_shadow.shad", ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *shadow_filename=new char[len];
    vsprintf( shadow_filename, "track%04d_shadow.shad", ap );
    
    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);   
    
    std::ofstream out;
    out.open( filename );
    
    delete[] filename;
    
    if( out == NULL )
    {
        return false;
    }
    
    // first line: what image file this RIB file should produce
    vasprintf( &filename, "track%04d.tiff", ap );
    
    char *shadow_filename;
    vasprintf( &shadow_filename, "track%04d_shadow.shad", ap );
    
    va_end(ap);
#endif
    // flatten
    std::vector<float> xs;
    for ( unsigned int i = 0; i < positions.size(); ++i )
    {
        xs.push_back( (float) positions[i][0] );
        xs.push_back( (float) positions[i][1] );
        xs.push_back( (float) positions[i][2] );
    }
    
    std::vector<Vec3f> normals;
    for ( unsigned int i = 0; i < positions.size(); ++i )
    {
        Vec3f n(0,0,0);
        for ( unsigned int j = 0; j < mesh.m_vertex_to_triangle_map[i].size(); ++j )
        {
            const Vec3st& tri = mesh.get_triangle( mesh.m_vertex_to_triangle_map[i][j] );
            Vec3d u = positions[ tri[1] ] - positions[ tri[0] ];
            Vec3d v = positions[ tri[2] ] - positions[ tri[0] ];
            Vec3d tn = normalized(cross(u, v));          
            n += Vec3f( (float)tn[0], (float)tn[1], (float)tn[2] ); 
        }
        normals.push_back( n / (float) mesh.m_vertex_to_triangle_map[i].size() );
    }
    
    
    out << "Display \"" << filename << "\" \"file\" \"rgb\"" << std::endl;
    delete[] filename;
    
    // next line: image format (width and height in pixels, pixel aspect ratio)
    out << "Format " << Gluvi::winwidth << " " << Gluvi::winheight << " 1" << std::endl;
    out << "PixelSamples 2 2" << std::endl;
    out << "Exposure 1.0 2.2" << std::endl;
    
    // then write out the camera specification
    Gluvi::camera->export_rib(out);
    
    // start the scene
    out << "WorldBegin\n";
    
    out << "LightSource \"ambientlight\" 1 \"intensity\" .3 \"lightcolor\" [1 1 1]\n";
    
    /*
     for ( unsigned int i = 0; i < lights.size(); ++i )
     {
     // compute location of light
     Vec3f light_pos( 0, 0, lights[i].dist );
     rotate( light_pos, lights[i].pitch, Vec3f(1,0,0) );
     rotate( light_pos, lights[i].heading, Vec3f(0,1,0) );
     light_pos += Vec3f( lights[i].target[0], lights[i].target[1], lights[i].target[2] );
     
     out << "LightSource \"singleshadowpoint\" 2 \"intensity\" 30 \"from\" [" << light_pos << "] \"shadowmap\" \"" << shadow_filename << "\"\n";
     }
     */
    
    //out << "LightSource \"distantlight\" 3 \"intensity\" 0.3 \"from\" [-5 -10 20] \"to\" [0 0 0]\n";
    
    out << "AttributeBegin\n";
    out << "  Color [0.6 0.6 0.2]\n";
    out << "  Opacity [1 1 1]\n";
    out << "  Surface \"matte\" \"Kd\" 1\n";
    
    const float plane_limit = 50.0f;
    const float plane_distance = 10.0f;
    char buf[256];
    sprintf( buf, "  Polygon \"P\" [-%f -%f -%f %f -%f -%f  %f %f -%f  -%f %f -%f]\n", 
            plane_limit, plane_limit, plane_distance, 
            plane_limit, plane_limit, plane_distance, 
            plane_limit, plane_limit, plane_distance, 
            plane_limit, plane_limit, plane_distance );
    out << buf;
    
    out << "AttributeEnd\n";
    
    out << "Color [0.3 0.3 0.9]\n";
    out << "Surface \"matte\" \n";
    
    write_ribfile( mesh, xs, out);
    
    out << " \"N\" [";
    for(unsigned int i=0; i<normals.size(); ++i)
    {
        out << normals[i] << "  ";
    }
    out << "]" << std::endl;
    
    // finish the scene
    out << "WorldEnd\n";
    
    out.flush();
    out.close();
    
    return true;
}

#endif

// ---------------------------------------------------------
//
// Write mesh in PBRT format
//
// ---------------------------------------------------------

bool write_pbrtfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, const char *filename_format, ...)
{

#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);
    std::ofstream output(filename, std::ofstream::binary);
    delete[] filename;
    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);
    std::ofstream output(filename, std::ofstream::binary);
    std::free(filename);
    va_end(ap);
#endif
    if(!output.good()) return false;
    return write_ribfile(mesh, x, output);
}

// ---------------------------------------------------------
//
// Write mesh in PBRT format
//
// ---------------------------------------------------------

bool write_pbrtfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, std::ostream &output)
{
    output<<"# generated by editmesh"<<std::endl;
    
    //output<<"\"integer nlevels\" [3]"<<std::endl;
    output<<"\"point P\" ["<<std::endl;
    for(unsigned int i=0; i<x.size(); ++i){
        output<<x[i]<<"  ";
        if(i%4==3 && i!=x.size()-1) output<<std::endl;
    }
    output<<"]"<<std::endl;
    output<<" \"integer indices\" ["<<std::endl;
    const std::vector<Vec3st>& tris = mesh.get_triangles();
    for(unsigned int i=0; i<tris.size(); ++i){
        output<<tris[i]<<"  ";
        if(i%6==5 && i!=tris.size()-1) output<<std::endl;
    }
    output<<"]"<<std::endl;
    
    return output.good();
}





