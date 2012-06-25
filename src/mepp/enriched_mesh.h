/*!
 * \file enriched_mesh.h
 * \brief Class: Enriched_mesh.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2012
 */

#ifndef _ENRICHED_MESH_
#define _ENRICHED_MESH_

#ifdef _MSC_VER
#pragma warning(disable: 4244 4267)
#endif

#ifdef _MSC_VER
// needed by OpenMesh to have M_PI etc defined
#if !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif
// some Mepp components want this on windows
#if !defined(NOMINMAX)
#define NOMINMAX
#endif
#endif

#include <iostream>

#include <QImage>
#include <QString>

#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#ifdef _MSC_VER
#include <windows.h>
#endif
#include <GL/glu.h>

// this needs to be done after all possible windows C header includes and before any Mepp components source includes
// (system C++ includes are supposed to be able to deal with this already):
// windows.h defines min and max macros which would make some Mepp components fail to compile.
#if defined(min) || defined(max)
#error The preprocessor symbols 'min' or 'max' are defined. If you are compiling on Windows, do #define NOMINMAX to prevent windows.h from defining these symbols.
#endif

/*#if defined(_MSC_VER)
#  undef min
#  undef max
#endif*/

template <typename M>
class Enriched_mesh
{
public:

  typedef M Mesh;

  // default constructor
  Enriched_mesh()
  {
  }
  
  // destructor
  ~Enriched_mesh() {}

  // open mesh
  /*virtual*/ bool open_mesh( const char* _filename, OpenMesh::IO::Options _opt );
  
  // load texture
  /*virtual*/ bool open_texture( const char *_filename );
  bool set_texture( QImage& _texsrc );

  Mesh& mesh() { return mesh_; }
  const Mesh& mesh() const { return mesh_; }

private:

  OpenMesh::IO::Options opt_;
  Mesh mesh_;
  
public: // TODO temp !!!
  OpenMesh::Vec3f bbMin, bbMax;
  float normal_scale_;
};

template <typename M>
bool 
Enriched_mesh<M>::open_mesh(const char* _filename, OpenMesh::IO::Options _opt)
{
  // load mesh
  // calculate normals
  
  mesh_.request_face_normals();
  mesh_.request_face_colors();
  mesh_.request_vertex_normals();
  mesh_.request_vertex_colors();
  mesh_.request_vertex_texcoords2D();
  
  std::cout << "--> Loading mesh from file '" << _filename << "'\n";
  if ( OpenMesh::IO::read_mesh(mesh_, _filename, _opt ))
  {
    // store read option
    opt_ = _opt;
    
    // update face and vertex normals     
    if ( ! opt_.check( OpenMesh::IO::Options::FaceNormal ) )
      mesh_.update_face_normals();
    else
      std::cout << "File provides face normals\n";
    
    if ( ! opt_.check( OpenMesh::IO::Options::VertexNormal ) )
      mesh_.update_vertex_normals();
    else
      std::cout << "File provides vertex normals\n";


    // check for possible color information
    if ( opt_.check( OpenMesh::IO::Options::VertexColor ) )
    {
      std::cout << "File provides vertex colors\n";
      //add_draw_mode("Colored Vertices");
    }
    else
      mesh_.release_vertex_colors();

    if ( _opt.check( OpenMesh::IO::Options::FaceColor ) )
    {
      std::cout << "File provides face colors\n";
      //add_draw_mode("Solid Colored Faces");
      //add_draw_mode("Smooth Colored Faces");
    }
    else
      mesh_.release_face_colors();

    if ( _opt.check( OpenMesh::IO::Options::VertexTexCoord ) )
      std::cout << "File provides texture coordinates\n";

    // bounding box
    typename Mesh::ConstVertexIter vIt(mesh_.vertices_begin());
    typename Mesh::ConstVertexIter vEnd(mesh_.vertices_end());

    bbMin = bbMax = OpenMesh::vector_cast<OpenMesh::Vec3f>(mesh_.point(vIt));
    
    for (size_t count=0; vIt!=vEnd; ++vIt, ++count)
    {
      bbMin.minimize( OpenMesh::vector_cast<OpenMesh::Vec3f>(mesh_.point(vIt)) );
      bbMax.maximize( OpenMesh::vector_cast<OpenMesh::Vec3f>(mesh_.point(vIt)) );
    }
        
    // set center and radius
    // TODO set_scene_pos( (bbMin+bbMax)*0.5, (bbMin-bbMax).norm()*0.5 );
    
    // for normal display
    normal_scale_ = (bbMax-bbMin).min()*0.05f;

    // loading done
    std::cout << "--> Loading done\n\n";
    return true;
  }
  std::cout << "--> Loading NOT done!\n\n";
  return false;
}

template <typename M>
bool Enriched_mesh<M>::open_texture( const char *_filename )
{
   QImage texsrc;
   QString fname = _filename;

   std::cout << "--> Loading texture from file '" << _filename << "'\n";
   if (texsrc.load( fname ))
   {
                std::cout << "--> Loading done\n\n";
                return true;
                //return set_texture( texsrc );
   }
   std::cout << "--> Loading NOT done!\n\n";
   return false;
}

#endif // _ENRICHED_MESH_
