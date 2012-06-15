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
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <ImporterCGAL.h>
#include <ExporterCGAL.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>


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
template <class HDS, class Mesh>
class Build_polyhedron : public CGAL::Modifier_base<HDS>
{
public:
    Build_polyhedron(Mesh& _mesh) : mesh_(_mesh) {}

    void operator()(HDS& hds)
	{
		// Postcondition: 'hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

        B.begin_surface(mesh_.n_vertices(), mesh_.n_faces(), mesh_.n_halfedges());

        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;

		// vertex
		typename Mesh::VertexIter vit;
		for(vit=mesh_.vertices_begin(); vit!=mesh_.vertices_end(); ++vit)
		{
			B.add_vertex(Point(mesh_.point(vit)[0], mesh_.point(vit)[1], mesh_.point(vit)[2]));
		}

		// facet
		typename Mesh::ConstFaceIter		fIt(mesh_.faces_begin()), fEnd(mesh_.faces_end());
		typename Mesh::ConstFaceVertexIter	fvIt;
		for (; fIt!=fEnd; ++fIt)
		{
			B.begin_facet();

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
		}

        B.end_surface();

		if (B.check_unconnected_vertices())
			B.remove_unconnected_vertices();
    }

private:
  Mesh& mesh_;
};

typedef CGAL::Simple_cartesian<float>		Kernel1;
typedef CGAL::Polyhedron_3<Kernel1>         Polyhedron1;
typedef Polyhedron1::HalfedgeDS             HalfedgeDS1;

typedef CGAL::Simple_cartesian<double>		Kernel2;
typedef CGAL::Polyhedron_3<Kernel2>         Polyhedron2;
typedef Polyhedron2::HalfedgeDS             HalfedgeDS2;

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
	OpenMesh::IO::ImporterCGAL/*<Mesh>*/ importer/*(_mesh)*/;
	std::cout << "--> Loading with ImporterCGAL\n";
	if ( OpenMesh::IO::IOManager().read("C:\\Users\\noname\\Desktop\\face.obj\\face.obj", importer, opt) )
		std::cout << "--> Loading done (vertices: " << importer.n_vertices() << " - faces: " << importer.n_faces() << ")\n\n";
	else
		std::cout << "--> Loading NOT done!\n\n";

	OpenMesh::IO::ExporterCGAL/*<Mesh>*/ exporter/*(_mesh)*/;
	std::cout << "--> Writing with ExporterCGAL\n";
	if ( OpenMesh::IO::IOManager().write("C:\\Users\\noname\\Desktop\\face_written.obj", exporter, opt) )
		std::cout << "--> Writing done (vertices: " << exporter.n_vertices() << " - faces: " << exporter.n_faces() << ")\n\n";
	else
		std::cout << "--> Writing NOT done!\n\n";
#endif

	std::cout << "--> OpenMeshToCGALConverter\n";
	Polyhedron1 P1;
    Build_polyhedron<HalfedgeDS1, MyMesh1> polyhedron1(mesh1);
    P1.delegate(polyhedron1);
	std::cout << "--> (vertices: " << P1.size_of_vertices() << " - faces: " << P1.size_of_facets() << " - edges: " << (P1.size_of_halfedges() >> 1) << ")\n\n";

	std::cout << "--> OpenMeshToCGALConverter\n";
	Polyhedron2 P2;
    Build_polyhedron<HalfedgeDS2, MyMesh2> polyhedron2(mesh2);
    P2.delegate(polyhedron2);
	std::cout << "--> (vertices: " << P2.size_of_vertices() << " - faces: " << P2.size_of_facets() << " - edges: " << (P2.size_of_halfedges() >> 1) << ")\n\n";

    return app.exec();
}