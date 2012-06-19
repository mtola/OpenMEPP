/*!
 * \file main.cpp
 * \brief Main file.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2012
 */

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
//#define NOMINMAX
#endif
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
//#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <ImporterCGAL.h>
#include <ExporterCGAL.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <enriched_polyhedron.h>

#include <OpenMesh/Core/IO/importer/ImporterT.hh>


#include <QApplication>
#include <QtGui>
//Q_IMPORT_PLUGIN(qjpeg)

//#include "mainwindow.hxx"

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


struct MyTraits : public OpenMesh::DefaultTraits
{
	HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
};
typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits>	MyMesh1;
//typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>	MyMesh1;

// Define my personal traits
struct MyTraits2 : OpenMesh::DefaultTraits
{
  // Let Point and Normal be a vector of doubles
  typedef OpenMesh::Vec3d Point;
  typedef OpenMesh::Vec3d Normal;

  // Already defined in OpenMesh::DefaultTraits
  // HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
  
  // Uncomment next line to disable attribute PrevHalfedge
  // HalfedgeAttributes( OpenMesh::Attributes::None );
  //
  // or
  //
  HalfedgeAttributes( 0 );
};
// Define my mesh with the new traits!
typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits2>	MyMesh2;
//typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits2>	MyMesh2;


template <typename M>
bool 
open_mesh( M& mesh_, const char* _filename, OpenMesh::IO::Options _opt )
{
  OpenMesh::IO::Options opt_; // mesh file options

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

	// TODO
	// bounding box

	// center and radius
	// normal scale for normal display / base point for displaying face normals

    // loading done
	std::cout << "--> Loading done\n\n";
    return true;
  }
  std::cout << "--> Loading NOT done!\n\n";
  return false;
}

template <typename M>
bool open_texture( const char *_filename )
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

/*!
 * \fn void loadStyleSheet()
 * \brief Load and apply StyleSheet (not used).
 */
void loadStyleSheet()
{
    // Let's use QFile and point to a resource...
    QFile data("./mepp.qss");
    QString style;

    // ...to open the file
    if (data.open(QFile::ReadOnly))
	{
        // QTextStream...
        QTextStream styleIn(&data);
        // ...read file to a string.
        style = styleIn.readAll();
        data.close();

        // We'll use qApp macro to get the QApplication pointer and set the style sheet application wide.
        qApp->setStyleSheet(style);
    }
}

// A modifier creating a polyhedron with the incremental builder.
template <class P, class HDS, class Mesh>
class Build_polyhedron : public CGAL::Modifier_base<HDS>
{
public:
    Build_polyhedron(Mesh& _mesh, P& _p) : mesh_(_mesh), p_(_p) {}

    void operator()(HDS& hds)
	{
		// Postcondition: 'hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

        B.begin_surface(mesh_.n_vertices(), mesh_.n_faces(), mesh_.n_halfedges());

        typedef typename HDS::Vertex			Vertex;
        typedef typename Vertex::Point			Point;

		typedef typename HDS::Vertex_handle		Vertex_handle;
		typedef typename HDS::Face_handle		Face_handle;

		// vertex
		typename Mesh::VertexIter vit;
		for(vit=mesh_.vertices_begin(); vit!=mesh_.vertices_end(); ++vit)
		{
			Vertex_handle vertex = B.add_vertex(Point(mesh_.point(vit)[0], mesh_.point(vit)[1], mesh_.point(vit)[2]));

			if (mesh_.has_vertex_colors())
			{
				OpenMesh::Vec3uc uc = mesh_.color(vit);
				vertex->color((float)uc[0]/255.0, (float)uc[1]/255.0, (float)uc[2]/255.0); // TODO unsigned char and alpha
			}
		}

		// facet
		typename Mesh::ConstFaceIter		fIt(mesh_.faces_begin()), fEnd(mesh_.faces_end());
		typename Mesh::ConstFaceVertexIter	fvIt;
		for (; fIt!=fEnd; ++fIt)
		{
			Face_handle face = B.begin_facet();

			// only for a triangle !!!
			/*fvIt = mesh_.cfv_iter(fIt.handle()); 
			B.add_vertex_to_facet(fvIt.handle().idx());
			++fvIt;
			B.add_vertex_to_facet(fvIt.handle().idx());
			++fvIt;
			B.add_vertex_to_facet(fvIt.handle().idx());*/

			for (fvIt=mesh_.cfv_iter(fIt.handle()); fvIt; ++fvIt)
				B.add_vertex_to_facet(fvIt.handle().idx());

			B.end_facet();

			if (mesh_.has_face_colors())
			{
				OpenMesh::Vec3uc uc = mesh_.color(fIt);
				face->color((float)uc[0]/255.0, (float)uc[1]/255.0, (float)uc[2]/255.0); // TODO unsigned char and alpha
			}
		}

        B.end_surface();

		if (B.check_unconnected_vertices())
			B.remove_unconnected_vertices();

		// ---

		p_.compute_normals();

		if (mesh_.has_vertex_colors())
			p_.has_vertex_colors(true);
		if (mesh_.has_face_colors())
			p_.has_face_colors(true);

		// TODO texture (VertexTexCoord)

		p_.compute_bounding_box();

		// TODO
		// center and radius
		// normal scale for normal display / base point for displaying face normals
    }

private:
  Mesh& mesh_;
  P& p_;
};

template <class Mesh, class P>
class Build_mesh
{
public:
    Build_mesh(P& _p) : p_(_p) {}

    void build(Mesh& mesh_)
	{
		typedef typename P::Vertex_const_iterator					VCI;
		typedef typename P::Facet_const_iterator					FCI;
		typedef typename P::Halfedge_around_facet_const_circulator	HFCC;

		OpenMesh::IO::ImporterT<Mesh> importerOpenMesh(mesh_);
		importerOpenMesh.reserve(p_.size_of_vertices(), 3*p_.size_of_vertices(), p_.size_of_facets());

		OpenMesh::Vec3f v;
		OpenMesh::VertexHandle vh;

		OpenMesh::IO::BaseImporter::VHandles vhandles;
		OpenMesh::FaceHandle fh;

		// output vertices
		for (VCI vi = p_.vertices_begin(); vi != p_.vertices_end(); ++vi)
		{
			v[0] = ::CGAL::to_double(vi->point().x());
			v[1] = ::CGAL::to_double(vi->point().y());
			v[2] = ::CGAL::to_double(vi->point().z());

			vh = importerOpenMesh.add_vertex(v);

			if (p_.has_vertex_colors())
			{
				mesh_.set_color(vh, OpenMesh::Vec3uc(128, 128, 128)); // TODO unsigned char and alpha
			}
		}

		// constructs an inverse index
		typedef CGAL::Inverse_index<VCI> Index;
		Index index(p_.vertices_begin(), p_.vertices_end());

		// output facets
		for (FCI fi = p_.facets_begin(); fi != p_.facets_end(); ++fi)
		{
			vhandles.clear(); // OpenMesh

			HFCC hc = fi->facet_begin();
			HFCC hc_end = hc;
			std::size_t n = circulator_size(hc);
			CGAL_assertion(n >= 3);

			// facet begin
			do
			{
				vhandles.push_back(OpenMesh::VertexHandle(index[VCI(hc->vertex())]));
				++hc;
			} while (hc != hc_end);
			// facet end

			fh = importerOpenMesh.add_face(vhandles);

			if (p_.has_face_colors())
			{
				mesh_.set_color(fh, OpenMesh::Vec3uc(128, 128, 128)); // TODO unsigned char and alpha
			}
		}
    }

private:
  P& p_;
};

typedef CGAL::Simple_cartesian<float>					Kernel1;
//typedef CGAL::Polyhedron_3<Kernel1>						Polyhedron1;
typedef Enriched_polyhedron<Kernel1, Enriched_items>	Polyhedron1;
typedef Polyhedron1::HalfedgeDS							HalfedgeDS1;

typedef CGAL::Simple_cartesian<double>					Kernel2;
//typedef CGAL::Polyhedron_3<Kernel2>						Polyhedron2;
typedef Enriched_polyhedron<Kernel2, Enriched_items>	Polyhedron2;
typedef Polyhedron2::HalfedgeDS							HalfedgeDS2;

/*!
 * \fn int main(int argc, char *argv[])
 * \brief Main function.
 *
 * \param argc .
 * \param argv .
 * \return 0 if normal.
 */
int main(int argc, char *argv[])
{
	QApplication app(argc, argv);

    //mainwindow window;
	QWidget window;

	loadStyleSheet();
	window.show();

	OpenMesh::IO::Options opt;

	MyMesh1 mesh1;	
	open_mesh<MyMesh1>(mesh1, "C:\\Users\\noname\\Desktop\\face.texture.ply\\face.ply", opt);
	std::cout << "--> (vertices: " << mesh1.n_vertices() << " - faces: " << mesh1.n_faces() << " - edges: " << mesh1.n_edges() << ")\n\n";

	MyMesh2 mesh2;
	open_mesh<MyMesh2>(mesh2, "C:\\Users\\noname\\Desktop\\face.obj\\face.obj", opt);
	std::cout << "--> (vertices: " << mesh2.n_vertices() << " - faces: " << mesh2.n_faces() << " - edges: " << mesh2.n_edges() << ")\n\n";

	open_texture<MyMesh1>("C:\\Users\\noname\\Desktop\\face.texture.ply\\face_skin_hi.bmp");
	open_texture<MyMesh2>("C:\\Users\\noname\\Desktop\\face.obj\\face_skin_hi.jpg");

#if(0)
	OpenMesh::IO::ImporterCGAL/*<Mesh>*/ importerCGAL/*(_mesh)*/;
	std::cout << "--> Loading with ImporterCGAL\n";
	if ( OpenMesh::IO::IOManager().read("C:\\Users\\noname\\Desktop\\face.obj\\face.obj", importerCGAL, opt) )
		std::cout << "--> Loading done (vertices: " << importerCGAL.n_vertices() << " - faces: " << importerCGAL.n_faces() << ")\n\n";
	else
		std::cout << "--> Loading NOT done!\n\n";

	OpenMesh::IO::ExporterCGAL/*<Mesh>*/ exporterCGAL/*(_mesh)*/;
	std::cout << "--> Writing with ExporterCGAL\n";
	if ( OpenMesh::IO::IOManager().write("C:\\Users\\noname\\Desktop\\face_written.obj", exporterCGAL, opt) )
		std::cout << "--> Writing done (vertices: " << exporterCGAL.n_vertices() << " - faces: " << exporterCGAL.n_faces() << ")\n\n";
	else
		std::cout << "--> Writing NOT done!\n\n";
#endif

	std::cout << "--> OpenMeshToCGALConverter\n";
	Polyhedron1 P1;
    Build_polyhedron<Polyhedron1, HalfedgeDS1, MyMesh1> build_polyhedron1(mesh1, P1);
    P1.delegate(build_polyhedron1);
	std::cout << "--> (vertices: " << P1.size_of_vertices() << " - faces: " << P1.size_of_facets() << " - edges: " << (P1.size_of_halfedges() >> 1) << ")\n\n";

	std::cout << "--> OpenMeshToCGALConverter\n";
	Polyhedron2 P2;
    Build_polyhedron<Polyhedron2, HalfedgeDS2, MyMesh2> build_polyhedron2(mesh2, P2);
    P2.delegate(build_polyhedron2);
	std::cout << "--> (vertices: " << P2.size_of_vertices() << " - faces: " << P2.size_of_facets() << " - edges: " << (P2.size_of_halfedges() >> 1) << ")\n\n";

	// ---

	std::cout << "--> CGALToOpenMeshConverter\n";
	MyMesh1 mesh1b;
	Build_mesh<MyMesh1, Polyhedron1> build_mesh1b(P1);
	build_mesh1b.build(mesh1b);
	std::cout << "--> (vertices: " << mesh1b.n_vertices() << " - faces: " << mesh1b.n_faces() << " - edges: " << mesh1b.n_edges() << ")\n\n";

	std::cout << "--> CGALToOpenMeshConverter\n";
	MyMesh2 mesh2b;
	Build_mesh<MyMesh2, Polyhedron2> build_mesh2b(P2);
	build_mesh2b.build(mesh2b);
	std::cout << "--> (vertices: " << mesh2b.n_vertices() << " - faces: " << mesh2b.n_faces() << " - edges: " << mesh2b.n_edges() << ")\n\n";

    return app.exec();
}