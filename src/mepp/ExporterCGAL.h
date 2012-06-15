/*!
 * \file ExporterCGAL.h
 * \brief Implements an exporter module for arbitrary CGAL meshes.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2012
 */
#ifndef __EXPORTERCGAL_H__
#define __EXPORTERCGAL_H__


//=== INCLUDES ================================================================
#include <OpenMesh/Core/IO/exporter/BaseExporter.hh>


//=== NAMESPACES ==============================================================
namespace OpenMesh {
namespace IO {


//=== IMPLEMENTATION ==========================================================

/**
 *  This class template provides an exporter module for CGAL meshes.
 */
//template <class Mesh>
class ExporterCGAL : public BaseExporter
{
public:

  // Constructor
  //ExporterCGAL(const Mesh& _mesh) : mesh_(_mesh) {}
  ExporterCGAL()
  {
	  vertices_ = faces_ = edges_ = 0;
  }

  // get vertex data
  Vec3f point(VertexHandle _vh) const 
  { 
    return Vec3f(0.0f, 0.0f, 0.0f); 
  }

  Vec3f normal(VertexHandle _vh) const 
  { 
    return Vec3f(0.0f, 0.0f, 0.0f);
  }

  Vec3uc color(VertexHandle _vh) const
  {
    return Vec3uc(0, 0, 0);
  }

  Vec4uc colorA(VertexHandle _vh) const
  {
    return Vec4uc(0, 0, 0, 0);
  }

  Vec2f texcoord(VertexHandle _vh) const
  {
	return Vec2f(0.0f, 0.0f);
  }
  
  // get edge data
  Vec3uc color(EdgeHandle _eh) const
  {
	return Vec3uc(0, 0, 0);
  }
  
  Vec4uc colorA(EdgeHandle _eh) const
  {
	return Vec4uc(0, 0, 0, 0);
  }

  // get face data
  unsigned int get_vhandles(FaceHandle _fh, std::vector<VertexHandle>& _vhandles) const
  {
    return 0;
  }

  Vec3f normal(FaceHandle _fh) const 
  { 
    return Vec3f(0.0f, 0.0f, 0.0f);
  }

  Vec3uc color(FaceHandle _fh) const 
  { 
    return Vec3uc(0, 0, 0);
  }

  Vec4uc colorA(FaceHandle _fh) const 
  { 
    return Vec4uc(0, 0, 0, 0);
  }

  // query number of faces, vertices, normals, texcoords
  size_t n_vertices()  const { return vertices_;/*mesh_.n_vertices();*/ }   
  size_t n_faces()     const { return faces_;/*mesh_.n_faces();*/ }
  size_t n_edges()     const { return edges_;/*mesh_.n_edges();*/ }

private:
  
	//const Mesh& mesh_;
	int vertices_;
	int faces_;
	int edges_;
};


//=============================================================================
} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif
//=============================================================================
