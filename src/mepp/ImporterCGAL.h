/*!
 * \file ImporterCGAL.h
 * \brief Implements an importer module for arbitrary CGAL meshes.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2012
 */
#ifndef __IMPORTERCGAL_H__
#define __IMPORTERCGAL_H__


//=== INCLUDES ================================================================
#include <OpenMesh/Core/IO/importer/BaseImporter.hh>


//=== NAMESPACES ==============================================================
namespace OpenMesh {
namespace IO {


//=== IMPLEMENTATION ==========================================================

/**
 *  This class template provides an importer module for CGAL meshes.
 */
//template <class Mesh>
class ImporterCGAL : public BaseImporter
{
public:

  typedef std::vector<VertexHandle>  VHandles;


  // Constructor
  //ImporterCGAL(Mesh& _mesh) : mesh_(_mesh) {}
  ImporterCGAL()
  {
	  vertices_ = faces_ = 0;
  }


  virtual VertexHandle add_vertex(const Vec3f& _point)
  {
	//std::cout << "---> add_vertex\n";
	return VertexHandle(vertices_++);
  }

  virtual FaceHandle add_face(const VHandles& _indices)
  {
	//std::cout << "---> add_face\n";
	return FaceHandle(faces_++);
  }

  // vertex attributes
  virtual void set_normal(VertexHandle _vh, const Vec3f& _normal)
  {
  }

  virtual void set_color(VertexHandle _vh, const Vec4uc& _color)
  {
  }

  virtual void set_color(VertexHandle _vh, const Vec3uc& _color)
  {
  }

  virtual void set_texcoord(VertexHandle _vh, const Vec2f& _texcoord)
  {
  }

  virtual void set_texcoord(HalfedgeHandle _heh, const Vec2f& _texcoord)
  {
  }

  // edge attributes  
  virtual void set_color(EdgeHandle _eh, const Vec4uc& _color)
  {
  }
  
  virtual void set_color(EdgeHandle _eh, const Vec3uc& _color)
  {
  }

  // face attributes
  virtual void set_normal(FaceHandle _fh, const Vec3f& _normal)
  {
  }

  virtual void set_color(FaceHandle _fh, const Vec3uc& _color)
  {
  }

  virtual void set_color(FaceHandle _fh, const Vec4uc& _color)
  {
  }

  virtual void add_face_texcoords( FaceHandle _fh, VertexHandle _vh, const std::vector<Vec2f>& _face_texcoords)
  {
  }

  virtual void set_face_texindex( FaceHandle _fh, int _texId )
  {
  }

  virtual void add_texture_information( int _id , std::string _name )
  {
  }

  // query number of faces, vertices, normals, texcoords
  size_t n_vertices()  const { return vertices_;/*mesh_.n_vertices();*/ }
  size_t n_faces()     const { return faces_;/*mesh_.n_faces();*/ }
  size_t n_edges()     const { return 0;/*mesh_.n_edges();*/ }

private:

	//Mesh& mesh_;
	int vertices_;
	int faces_;
};

//=============================================================================
} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif
