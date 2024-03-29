/*!
 * \file enriched_polyhedron.h
 * \brief Class: Enriched_polyhedron.
 * \author Le-Jeng Shiue, Pierre Alliez, Radu Ursu, Lutz Kettner and Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2004-2013
 */

#ifndef _ENRICHED_POLYHEDRON_
#define _ENRICHED_POLYHEDRON_

// CGAL stuff
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include <fstream>
#include <list>

#ifdef _MSC_VER
#include <windows.h>
#endif
#include <GL/glu.h>

// compute facet normal 
struct Facet_normal // (functor)
{
  template <class Facet>
  void operator()(Facet& f)
  {
    typename Facet::Normal_3 sum = CGAL::NULL_VECTOR;
    typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
    do
    {
      typename Facet::Normal_3 normal = CGAL::cross_product(
        h->next()->vertex()->point() - h->vertex()->point(),
        h->next()->next()->vertex()->point() - h->next()->vertex()->point());
      double sqnorm = CGAL::to_double(normal * normal);
      if(sqnorm != 0)
        normal = normal / (float)std::sqrt(sqnorm);
      sum = sum + normal;
    }
    while(++h != f.facet_begin());
    double sqnorm = CGAL::to_double(sum * sum);
    if(sqnorm != 0.0)
      f.normal() = sum / std::sqrt(sqnorm);
    else
    {
      f.normal() = CGAL::NULL_VECTOR;
      //TRACE("degenerate face\n");
    }
  }
};

// compute vertex normal 
struct Vertex_normal // (functor)
{
    template <class Vertex>
    void operator()(Vertex& v)
    {
        typename Vertex::Normal_3 normal = CGAL::NULL_VECTOR;
        typename Vertex::Halfedge_around_vertex_const_circulator pHalfedge = v.vertex_begin();
        typename Vertex::Halfedge_around_vertex_const_circulator begin = pHalfedge;
        CGAL_For_all(pHalfedge,begin) 
          if(!pHalfedge->is_border())
            normal = normal + pHalfedge->facet()->normal();
        double sqnorm = CGAL::to_double(normal * normal);
        if(sqnorm != 0.0f)
          v.normal() = normal / (float)std::sqrt(sqnorm);
        else
          v.normal() = CGAL::NULL_VECTOR;
    }
};

// a refined facet with a normal and a tag
template <class Refs, class T, class P, class Norm>
class Enriched_facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
private:

  // tag
  int m_tag; 

  // normal
  Norm m_normal;

  // color
  unsigned char m_color[3];

public:

  // life cycle
  // no constructors to repeat, since only
  // default constructor mandatory
  Enriched_facet()
  {
	color(127/*0.5f*/, 127/*0.5f*/, 127/*0.5f*/);
  }

  // tag
  const int& tag() { return m_tag; }
  void tag(const int& t)  { m_tag = t; }

  // normal
  typedef Norm Normal_3;
  Normal_3& normal() { return m_normal; }
  const Normal_3& normal() const { return m_normal; }

  // color
  unsigned char color(int index) { return m_color[index]; }
  void color(unsigned char r, unsigned char g, unsigned char b) { m_color[0] = r; m_color[1] = g; m_color[2] = b; }
};

// a refined halfedge with a general tag and 
// a binary tag to indicate wether it belongs 
// to the control mesh or not
template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Enriched_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
private:

  // tag
  int m_tag; 

  // option for edge superimposing
  //bool m_control_edge; 

public:

  // life cycle
  Enriched_halfedge()
  {
    //m_control_edge = true;
  }

  // tag
  const int& tag() const { return m_tag;  }
  int& tag() { return m_tag;  }
  void tag(const int& t)  { m_tag = t; }

  // control edge 
  //bool& control_edge()  { return m_control_edge; }
  //const bool& control_edge()  const { return m_control_edge; }
  //void control_edge(const bool& flag) { m_control_edge  = flag; }
};

// a refined vertex with a normal and a tag
template <class Refs, class T, class P, class Norm>
class Enriched_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
private:

  // tag
  int m_tag; 

  // normal
  Norm m_normal;

  // color
  unsigned char m_color[3];

  // texture
  float m_texture_coordinates[2];

public:
  // life cycle
  Enriched_vertex()
  {
	color(127/*0.5f*/, 127/*0.5f*/, 127/*0.5f*/);
	texcoord2D(0.f, 0.f);
  }

  // repeat mandatory constructors
  Enriched_vertex(const P& pt)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt)
  {
  }

  // normal
  typedef Norm Normal_3;
  Normal_3& normal() { return m_normal; }
  const Normal_3& normal() const { return m_normal; }

  // tag
  int& tag() {  return m_tag; }
  const int& tag() const {  return m_tag; }
  void tag(const int& t)  { m_tag = t; }

  // color
  unsigned char color(int index) { return m_color[index]; }
  void color(unsigned char r, unsigned char g, unsigned char b) { m_color[0] = r; m_color[1] = g; m_color[2] = b; }

  // texture
  float texcoord2D(int index) { return m_texture_coordinates[index]; }
  void texcoord2D(float u, float v) { m_texture_coordinates[0] = u; m_texture_coordinates[1] = v; }
};

// A redefined items class for the Polyhedron_3 
// with a refined vertex class that contains a 
// member for the normal vector and a refined
// facet with a normal vector instead of the 
// plane equation (this is an alternative 
// solution instead of using 
// Polyhedron_traits_with_normals_3).
struct Enriched_items : public CGAL::Polyhedron_items_3
{
    // wrap vertex
    template <class Refs, class Traits>
    struct Vertex_wrapper // is a vertex with a reference to an incident halfedge and it stores a point of type Point
    {
        typedef typename Traits::Point_3  Point;
        typedef typename Traits::Vector_3 Normal;
        typedef Enriched_vertex<Refs,
                          CGAL::Tag_true,
                          Point,
                          Normal> Vertex;
    };

    // wrap face
    template <class Refs, class Traits>
    struct Face_wrapper // is a face with a reference to an incident halfedge, it stores a normal vector instead of a plane equation of type Plane
    {
        typedef typename Traits::Point_3  Point;
        typedef typename Traits::Vector_3 Normal;
        typedef Enriched_facet<Refs,
                         CGAL::Tag_true,
                         Point,
                         Normal> Face;
    };

    // wrap halfedge
    template <class Refs, class Traits> // in all cases, a reference to the next halfedge and to the opposite halfedge is supported
    struct Halfedge_wrapper // here, a reference to the previous halfedge, to the incident vertex and to the incident face is also supported
    {
        typedef typename Traits::Vector_3 Normal;
        typedef Enriched_halfedge<Refs,
                            CGAL::Tag_true,
                            CGAL::Tag_true,
                            CGAL::Tag_true,
                            Normal> Halfedge;
    };
};

//*********************************************************
template <class kernel, class items>
class Enriched_polyhedron : public CGAL::Polyhedron_3<kernel,items>
{
public :
  typedef typename kernel::FT FT;
  typedef typename kernel::Point_3 Point;
  typedef typename kernel::Vector_3 Vector;
  typedef typename kernel::Iso_cuboid_3 Iso_cuboid;
  
  // ---
  
  typedef typename Enriched_polyhedron::Facet_handle Facet_handle;
  typedef typename Enriched_polyhedron::Vertex_handle Vertex_handle;
  typedef typename Enriched_polyhedron::Halfedge_handle Halfedge_handle;
  
  typedef typename Enriched_polyhedron::Facet_iterator Facet_iterator;
  typedef typename Enriched_polyhedron::Vertex_iterator Vertex_iterator;
  typedef typename Enriched_polyhedron::Halfedge_iterator Halfedge_iterator;
  typedef typename Enriched_polyhedron::Edge_iterator Edge_iterator;
  typedef typename Enriched_polyhedron::Point_iterator Point_iterator;
  
  typedef typename Enriched_polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
  typedef typename Enriched_polyhedron::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
  
  typedef typename Enriched_polyhedron::HalfedgeDS HalfedgeDS;
  typedef typename HalfedgeDS::Face Facet;
  typedef typename Facet::Normal_3 Normal;

private :

  // bounding box
  Iso_cuboid m_bbox;

  // type
  bool m_pure_quad;
  bool m_pure_triangle;

  // color
  bool m_vertex_color;
  bool m_face_color;

  // texture
  bool m_vertex_texcoords2D;

private:

  GLuint tex_id_;

  Vector bbMin_, bbMax_;
  float normal_scale_;

  Vector center_;
  float radius_;

public :

  // life cycle
  Enriched_polyhedron() 
  {
    m_pure_quad = false;
    m_pure_triangle = false;

	// color
	m_vertex_color = false;
    m_face_color = false;

	// texture
	m_vertex_texcoords2D = false;

	normal_scale_ = radius_ = 0.;
  }
  virtual ~Enriched_polyhedron() 
  {
  }

  // type
  bool is_pure_triangle() { return m_pure_triangle; }
  bool is_pure_quad() { return m_pure_quad; }

  // color
  bool has_vertex_colors() { return m_vertex_color; }
  bool has_face_colors() { return m_face_color; }

  void has_vertex_colors(const bool& flag) { m_vertex_color = flag; }
  void has_face_colors(const bool& flag) { m_face_color = flag; }

  // texture
  bool has_vertex_texcoords2D() { return m_vertex_texcoords2D; }
  void has_vertex_texcoords2D(const bool& flag) { m_vertex_texcoords2D = flag; }

  GLuint& tex_id() { return tex_id_; }
  void tex_id(GLuint t) { tex_id_ = t; }

  // ---
  Vector& bbMin() { return bbMin_; }
  void bbMin(OpenMesh::Vec3f min) { bbMin_ = min; }

  Vector& bbMax() { return bbMax_; }
  void bbMax(OpenMesh::Vec3f max) { bbMax_ = max; }

  float& normal_scale() { return normal_scale_; }
  void normal_scale(float n) { normal_scale_ = n; }

  void center(Vector c) { center_ = c; }

  float& radius() { return radius_; }
  void radius(float r) { radius_ = r; }

  // normals (per facet, then per vertex)
  void compute_normals_per_facet()
  {
    std::for_each(this->facets_begin(),this->facets_end(),Facet_normal());
  }
  void compute_normals_per_vertex()
  {
    std::for_each(this->vertices_begin(),this->vertices_end(),Vertex_normal());
  }
  void compute_normals()
  {
    compute_normals_per_facet();
    compute_normals_per_vertex();
  }

  // bounding box
  Iso_cuboid& bbox() { return m_bbox; }
  const Iso_cuboid bbox() const { return m_bbox; }

  // compute bounding box
  void compute_bounding_box()
  {
    if(this->size_of_vertices() == 0)
    {
      /*ASSERT*/CGAL_assertion(false);
      return;
    }

    FT xmin,xmax,ymin,ymax,zmin,zmax;
    Vertex_iterator pVertex = this->vertices_begin();
    xmin = xmax = pVertex->point().x();
    ymin = ymax = pVertex->point().y();
    zmin = zmax = pVertex->point().z();
    for(;
        pVertex !=  this->vertices_end();
        pVertex++)
    {
      const Point& p = pVertex->point();

      xmin =  std::min(xmin,p.x());
      ymin =  std::min(ymin,p.y());
      zmin =  std::min(zmin,p.z());

      xmax =  std::max(xmax,p.x());
      ymax =  std::max(ymax,p.y());
      zmax =  std::max(zmax,p.z());
    }
    m_bbox = Iso_cuboid(xmin,ymin,zmin,
                        xmax,ymax,zmax);

	bbMin_ = Vector(xmin,ymin,zmin);
	bbMax_ = Vector(xmax,ymax,zmax);
  }

  // bounding box
  FT xmin() { return m_bbox.xmin(); }
  FT xmax() { return m_bbox.xmax(); }
  FT ymin() { return m_bbox.ymin(); }
  FT ymax() { return m_bbox.ymax(); }
  FT zmin() { return m_bbox.zmin(); }
  FT zmax() { return m_bbox.zmax(); }

  // copy bounding box
  void copy_bounding_box(Enriched_polyhedron<kernel,items> *pMesh)
  {
    m_bbox = pMesh->bbox();
  }

  // degree of a face
  static unsigned int degree(Facet_handle pFace)
  {
    return CGAL::circulator_size(pFace->facet_begin());    
  }

  // valence of a vertex
  static unsigned int valence(Vertex_handle pVertex)
  {
    return CGAL::circulator_size(pVertex->vertex_begin());
  }

  // check wether a vertex is on a boundary or not
  static bool is_border(Vertex_handle pVertex)
  {
    Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
    if(pHalfEdge == NULL) // isolated vertex
      return true;
    Halfedge_around_vertex_circulator d = pHalfEdge;
    CGAL_For_all(pHalfEdge,d)
      if(pHalfEdge->is_border())
        return true;
    return false;
  }

  // get any border halfedge attached to a vertex
  Halfedge_handle get_border_halfedge(Vertex_handle pVertex)
  {
    Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator d = pHalfEdge;
    CGAL_For_all(pHalfEdge,d)
      if(pHalfEdge->is_border())
        return pHalfEdge;
    return NULL;
  }

  // tag all halfedges
  void tag_halfedges(const int tag)
  {
    for(Halfedge_iterator pHalfedge = this->halfedges_begin();
        pHalfedge != this->halfedges_end();
        pHalfedge++)
      pHalfedge->tag(tag);
  }

  // tag all facets
  void tag_facets(const int tag)
  {
    for(Facet_iterator pFacet = this->facets_begin();
        pFacet  != this->facets_end();
        pFacet++)
      pFacet->tag(tag);
  }

  // set index for all vertices
  void set_index_vertices()
  {
    int index = 0;
    for(Vertex_iterator pVertex = this->vertices_begin();
        pVertex != this->vertices_end();
        pVertex++)
      pVertex->tag(index++);
  }

  // is pure degree ?
  bool is_pure_degree(unsigned int d)
  {
    for(Facet_iterator pFace  = this->facets_begin();
        pFace != this->facets_end();
        pFace++)
      if(degree(pFace) != d)
        return false;
    return true;
  }

  // compute type
  void compute_type()
  {
    m_pure_quad = is_pure_degree(4);
    m_pure_triangle = is_pure_degree(3);
  }

  // compute facet center
  void compute_facet_center(Facet_handle pFace,
                            Point& center)
  {
    Halfedge_around_facet_circulator pHalfEdge = pFace->facet_begin();
    Halfedge_around_facet_circulator end = pHalfEdge;
    Vector vec(0.0,0.0,0.0);
    int degree = 0;
    CGAL_For_all(pHalfEdge,end)
    {
      vec = vec + (pHalfEdge->vertex()->point()-CGAL::ORIGIN);
      degree++;
    }
    center = CGAL::ORIGIN + (vec/(/*kernel::*/FT)degree);
  }

  // compute average edge length around a vertex
  FT average_edge_length_around(Vertex_handle pVertex)
  {
    FT sum = 0.0;
    Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator end = pHalfEdge;
    Vector vec(0.0,0.0,0.0);
    int degree = 0;
    CGAL_For_all(pHalfEdge,end)
    {
      Vector vec = pHalfEdge->vertex()->point()-
                   pHalfEdge->opposite()->vertex()->point();
      sum += std::sqrt(CGAL::to_double(vec*vec));
      degree++;
    }
    return sum / (FT) degree;
  }

  // draw using OpenGL commands (display lists)
  void gl_draw(bool smooth_shading,
               bool use_normals)
  {
    // draw polygons
    Facet_iterator pFacet = this->facets_begin();
    for(;pFacet != this->facets_end();pFacet++)
    {
      // begin polygon assembly
      ::glBegin(GL_POLYGON);
        gl_draw_facet(pFacet,smooth_shading,use_normals);
      ::glEnd(); // end polygon assembly
    }
    glFlush();
  }

  void gl_draw_facet(Facet_handle pFacet,
                      bool smooth_shading,
                      bool use_normals)
  {
    // one normal per face
    if(use_normals && !smooth_shading)
    {
      const /*Facet::*/Normal/*_3*/& normal = pFacet->normal();
      ::glNormal3d(CGAL::to_double(normal[0]),CGAL::to_double(normal[1]),CGAL::to_double(normal[2]));
    }

    // revolve around current face to get vertices
    Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
    do
    {
      // one normal per vertex
      if(use_normals && smooth_shading)
      {
        const /*Facet::*/Normal/*_3*/& normal = pHalfedge->vertex()->normal();
        ::glNormal3d(CGAL::to_double(normal[0]),CGAL::to_double(normal[1]),CGAL::to_double(normal[2]));
      }

      // polygon assembly is performed per vertex
      const Point& point  = pHalfedge->vertex()->point();
      ::glVertex3d(CGAL::to_double(point[0]),CGAL::to_double(point[1]),CGAL::to_double(point[2]));
    }
    while(++pHalfedge != pFacet->facet_begin());
  }

  // superimpose edges
  void superimpose_edges(bool skip_ordinary_edges = true,
                         bool skip_control_edges = false,
                         bool voronoi_edge = false)
  {
    ::glBegin(GL_LINES);
    for(Edge_iterator h = this->edges_begin();
        h != this->edges_end();
        h++)
    {
      // ignore this edges
      if(skip_ordinary_edges && !h->control_edge())
        continue;

      // ignore control edges
      if(skip_control_edges && h->control_edge())
        continue;

      if(voronoi_edge)
      {
        Facet_handle pFace1 = h->facet();
        Facet_handle pFace2 = h->opposite()->facet();
        if(pFace1 == NULL || pFace2 == NULL)
          continue;

        const Point &p1 = h->vertex()->point();
        const Point &p2 = h->next()->vertex()->point();
        const Point &p3 = h->next()->next()->vertex()->point();

        kernel k;
        Point d1 = k.construct_circumcenter_3_object()(p1,p2,p3);
        ::glVertex3d(CGAL::to_double(d1[0]),CGAL::to_double(d1[1]),CGAL::to_double(d1[2]));

        const Point &pp1 = h->opposite()->vertex()->point();
        const Point &pp2 = h->opposite()->next()->vertex()->point();
        const Point &pp3 = h->opposite()->next()->next()->vertex()->point();
        Point d2 = k.construct_circumcenter_3_object()(pp1,pp2,pp3);
        ::glVertex3d(CGAL::to_double(d2[0]),CGAL::to_double(d2[1]),CGAL::to_double(d2[2]));
      }
      else
      {
        // assembly and draw line segment
        const Point& p1 = h->prev()->vertex()->point();
        const Point& p2 = h->vertex()->point();
        ::glVertex3d(CGAL::to_double(p1[0]),CGAL::to_double(p1[1]),CGAL::to_double(p1[2]));
        ::glVertex3d(CGAL::to_double(p2[0]),CGAL::to_double(p2[1]),CGAL::to_double(p2[2]));
      }
    }
    ::glEnd();
  }

  // superimpose vertices
  void superimpose_vertices()
  {
    ::glBegin(GL_POINTS);
    for(Point_iterator pPoint = this->points_begin();
        pPoint != this->points_end();
        pPoint++)
      ::glVertex3d(CGAL::to_double(pPoint->x()),CGAL::to_double(pPoint->y()),CGAL::to_double(pPoint->z()));
    ::glEnd(); // // end point assembly
  }

  // superimpose vertices
  void superimpose_spheres(double scale)
  {
    GLUquadricObj* pQuadric = gluNewQuadric();
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  
    for(Vertex_iterator pVertex = this->vertices_begin();
        pVertex !=  this->vertices_end();
        pVertex++)
    {
      ::glPushMatrix();
      double radius = CGAL::to_double(average_edge_length_around(pVertex));
      ::glTranslated(CGAL::to_double(pVertex->point().x()),
                     CGAL::to_double(pVertex->point().y()),
                     CGAL::to_double(pVertex->point().z()));
      ::gluSphere(pQuadric,scale*radius,24,24); 
      ::glPopMatrix();
    }
    gluDeleteQuadric(pQuadric);
  }

  // write in obj file format (OBJ).
  void write_obj(char *pFilename,
                  int incr  = 1) // 1-based by default
  {
    std::ofstream stream(pFilename);

    // output vertices
    for(Point_iterator pPoint = this->points_begin();
        pPoint != this->points_end(); 
        pPoint++) 
      stream << 'v' << ' ' << pPoint->x() << ' ' <<
                              pPoint->y() << ' ' <<
                              pPoint->z() << std::endl;		// CGAL::to_double here ???

    // precompute vertex indices
    this->set_index_vertices(); 

    // output facets
    for(Facet_iterator pFacet = this->facets_begin();
        pFacet != this->facets_end(); 
        pFacet++) 
    {
      stream << 'f';
      Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
      do 
        stream << ' ' << pHalfedge->vertex()->tag()+incr;
      while(++pHalfedge != pFacet->facet_begin());
      stream << std::endl;
    }  
  }

  // draw bounding box
  void gl_draw_bounding_box()
  {
    ::glBegin(GL_LINES);

	// along x axis
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));

	// along y axis
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));

	// along z axis
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
	::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));

    ::glEnd();
  }

  // count #boundaries
  unsigned int nb_boundaries()
  {
    unsigned int nb = 0;
    tag_halfedges(0);
    for(Halfedge_iterator he = this->halfedges_begin();
        he != this->halfedges_end();
        he++)
    {
      if(he->is_border() && he->tag() == 0)
      {
        nb++;
        Halfedge_handle curr = he;
        do
        {
          curr  = curr->next();
          curr->tag(1);
        }
        while(curr != he);
      }
    }
    return nb;
  }

  // tag component 
  void tag_component(Facet_handle pSeedFacet,
                     const int tag_free,
                     const int tag_done)
  {
    pSeedFacet->tag(tag_done);
    std::list<Facet_handle> facets;
    facets.push_front(pSeedFacet);
    while(!facets.empty())
    {
      Facet_handle pFacet = facets.front();
      facets.pop_front();
      pFacet->tag(tag_done);
      Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
      Halfedge_around_facet_circulator end = pHalfedge;
      CGAL_For_all(pHalfedge,end)
      {
        Facet_handle pNFacet = pHalfedge->opposite()->facet();
        if(pNFacet != NULL && pNFacet->tag() == tag_free)
        {
          facets.push_front(pNFacet);
          pNFacet->tag(tag_done);
        }
      }
    }
  }

  // count #components
  unsigned int nb_components()
  {
    unsigned int nb = 0;
    tag_facets(0);
    for(Facet_iterator pFacet = this->facets_begin();
        pFacet != this->facets_end();
        pFacet++)
    {
      if(pFacet->tag() == 0)
      {
        nb++;
        tag_component(pFacet,0,1);
      }
    }
    return nb;
  }

  // compute the genus
  // V - E + F + B = 2 (C - G)
  // C -> #connected components
  // G : genus
  // B : #boundaries
  int genus()
  {
    int c = nb_components();
    int b = nb_boundaries();
    int v = this->size_of_vertices();
    int e = this->size_of_halfedges()/2;
    int f = this->size_of_facets();
    return genus(c,v,f,e,b);
  }
  int genus(int c,
            int v,
            int f,
            int e,
            int b)
  {
    return (2*c+e-b-f-v)/2;
  }

  // ---

	void triangulate()
	{
		Facet_iterator f = this->facets_begin();
		Facet_iterator f2 = this->facets_begin();
		do //for (; f != this->facets_end(); f++)
		{
			f = f2;
			if (f == this->facets_end())
			{
				break;
			}
			f2++;

			if (!(f->is_triangle()))
			{
				int num = (int)(f->facet_degree() - 3);
				Halfedge_handle h = f->halfedge();

				h = this->make_hole(h);

				Halfedge_handle g = h->next();
				g = g->next();

				g = this->add_facet_to_border (h, g);

				num--;
				while (num != 0)
				{
					g = g->opposite();
					g = g->next();
					g = this->add_facet_to_border (h, g);
					num--;
				}

				this->fill_hole(h);
			}

		} while (true);

		this->compute_normals();
		this->compute_type();
	}
};

#endif // _ENRICHED_POLYHEDRON_
