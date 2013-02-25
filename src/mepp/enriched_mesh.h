/*!
 * \file enriched_mesh.h
 * \brief Class: Enriched_mesh.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2013
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

#include <QGLWidget>

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
	  normal_scale_ = radius_ = 0.;
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

  // texture
  GLuint& tex_id() { return tex_id_; }
  void tex_id(GLuint t) { tex_id_ = t; }

  // ---
  OpenMesh::Vec3f& bbMin() { return bbMin_; }
  void bbMin(OpenMesh::Vec3f min) { bbMin_ = min; }

  OpenMesh::Vec3f& bbMax() { return bbMax_; }
  void bbMax(OpenMesh::Vec3f max) { bbMax_ = max; }

  float& normal_scale() { return normal_scale_; }
  void normal_scale(float n) { normal_scale_ = n; }

  void center(OpenMesh::Vec3f c) { center_ = c; }

  float& radius() { return radius_; }
  void radius(float r) { radius_ = r; }

private:

  OpenMesh::IO::Options opt_;
  Mesh mesh_;

private:

  GLuint tex_id_;

  OpenMesh::Vec3f bbMin_, bbMax_;
  float normal_scale_;

  OpenMesh::Vec3f center_;
  float radius_;
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

    if ( opt_.check( OpenMesh::IO::Options::FaceColor ) )
    {
      std::cout << "File provides face colors\n";
      //add_draw_mode("Solid Colored Faces");
      //add_draw_mode("Smooth Colored Faces");
    }
    else
      mesh_.release_face_colors();

    if ( opt_.check( OpenMesh::IO::Options::VertexTexCoord ) )
      std::cout << "File provides texture coordinates\n";
	else
      mesh_.release_vertex_texcoords2D();

    // bounding box
    typename Mesh::ConstVertexIter vIt(mesh_.vertices_begin());
    typename Mesh::ConstVertexIter vEnd(mesh_.vertices_end());

    bbMin_ = bbMax_ = OpenMesh::vector_cast<OpenMesh::Vec3f>(mesh_.point(vIt));
    
    for (size_t count=0; vIt!=vEnd; ++vIt, ++count)
    {
      bbMin_.minimize( OpenMesh::vector_cast<OpenMesh::Vec3f>(mesh_.point(vIt)) );
      bbMax_.maximize( OpenMesh::vector_cast<OpenMesh::Vec3f>(mesh_.point(vIt)) );
    }
        
    // set center and radius ( set_scene_pos( , ) )
    center_ = (bbMin_+bbMax_)*0.5;
    radius_ = (bbMin_-bbMax_).norm()*0.5;
    
    // for normal display
    normal_scale_ = (bbMax_-bbMin_).min()*0.05f;

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
	if (set_texture( texsrc ))
	{
		std::cout << "--> Texture loading done\n\n";
		return true;
	}

   std::cout << "--> Texture loading NOT done!\n\n";
   return false;
}

template <typename M>
bool Enriched_mesh<M>::set_texture( QImage& _texsrc )
{
  if ( !opt_.vertex_has_texcoord() )
    return false;
   
  {
    // adjust texture size: 2^k * 2^l
    int tex_w, w( _texsrc.width()  );
    int tex_h, h( _texsrc.height() );

    for (tex_w=1; tex_w <= w; tex_w <<= 1) {};
    for (tex_h=1; tex_h <= h; tex_h <<= 1) {};
    tex_w >>= 1;
    tex_h >>= 1;
    _texsrc = _texsrc.scaled( tex_w, tex_h, Qt::IgnoreAspectRatio, Qt::SmoothTransformation );
  }

  QImage texture( QGLWidget::convertToGLFormat ( _texsrc ) );
  
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_UNPACK_SKIP_ROWS,   0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
  glPixelStorei(GL_UNPACK_ALIGNMENT,   1);
  glPixelStorei(GL_PACK_ROW_LENGTH,    0);
  glPixelStorei(GL_PACK_SKIP_ROWS,     0);
  glPixelStorei(GL_PACK_SKIP_PIXELS,   0);
  glPixelStorei(GL_PACK_ALIGNMENT,     1);    
  
  if ( tex_id_ > 0 )
  {
    glDeleteTextures(1, &tex_id_);
  }
  glGenTextures(1, &tex_id_);
  glBindTexture(GL_TEXTURE_2D, tex_id_);
    
  // glTexGenfv( GL_S, GL_SPHERE_MAP, 0 );
  // glTexGenfv( GL_T, GL_SPHERE_MAP, 0 );
    
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);      
  
  glTexImage2D(GL_TEXTURE_2D,       // target
	       0,                   // level
	       GL_RGBA,             // internal format
	       texture.width(),     // width  (2^n)
	       texture.height(),    // height (2^m)
	       0,                   // border
	       GL_RGBA,             // format
	       GL_UNSIGNED_BYTE,    // type
	       texture.bits() );    // pointer to pixels

  return true;
}

#endif // _ENRICHED_MESH_
