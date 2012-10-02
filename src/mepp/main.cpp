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
#define NOMINMAX
#endif
#endif

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
//#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <enriched_mesh.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <enriched_polyhedron.h>

#include <mesh_to_polyhedron.h>

//#include <ImporterCGAL.h>
//#include <ExporterCGAL.h>

#include <QApplication>
#include <QtGui>
//Q_IMPORT_PLUGIN(qjpeg)

//#include "mainwindow.hxx"


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


struct MyTraits : public OpenMesh::DefaultTraits
{
        HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
};
typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits>       MyMesh1;
//typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>      MyMesh1;

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
typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits2>      MyMesh2;
//typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits2>     MyMesh2;


typedef CGAL::Simple_cartesian<float>					Kernel1;
//typedef CGAL::Polyhedron_3<Kernel1>						Polyhedron1;
typedef Enriched_polyhedron<Kernel1, Enriched_items>	Polyhedron1;
typedef Polyhedron1::HalfedgeDS							HalfedgeDS1;

typedef CGAL::Simple_cartesian<double>					Kernel2;
//typedef CGAL::Polyhedron_3<Kernel2>						Polyhedron2;
typedef Enriched_polyhedron<Kernel2, Enriched_items>	Polyhedron2;
typedef Polyhedron2::HalfedgeDS							HalfedgeDS2;


typedef Enriched_mesh<MyMesh1> Mesh1;
typedef Enriched_mesh<MyMesh2> Mesh2;


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

	OpenMesh::IO::Options opt1; opt1 += OpenMesh::IO::Options::VertexTexCoord; opt1 += OpenMesh::IO::Options::VertexColor; opt1 += OpenMesh::IO::Options::FaceColor;
	OpenMesh::IO::Options opt2; //opt2 += OpenMesh::IO::Options::VertexColor;

	Mesh1 M1;
    //M1.open_mesh("C:\\Users\\noname\\Desktop\\face.texture.ply\\face.ply", opt1);
	//M1.open_mesh("/home/mepp/Desktop/face.texture.ply/face.ply", opt1);
	//M1.open_mesh("mesh_c.off", opt1);
	M1.open_mesh("C:\\Users\\noname\\Desktop\\_3dvia obj mepp_\\kip.obj\\kip_vt.ply", opt1);
	//M1.open_mesh("C:\\Users\\noname\\Desktop\\meshes\\rgb_monkey.obj", opt1);
	std::cout << "--> (vertices: " << M1.mesh().n_vertices() << " - faces: " << M1.mesh().n_faces() << " - edges: " << M1.mesh().n_edges() << ")\n\n";

	Mesh2 M2;
    M2.open_mesh("C:\\Users\\noname\\Desktop\\face.obj\\face.obj", opt2);
    //M2.open_mesh("/home/mepp/Desktop/face.obj/face.obj", opt2);
	//M2.open_mesh("mesh_t.off", opt2);
	std::cout << "--> (vertices: " << M2.mesh().n_vertices() << " - faces: " << M2.mesh().n_faces() << " - edges: " << M2.mesh().n_edges() << ")\n\n";

	M1.open_texture("C:\\Users\\noname\\Desktop\\face.texture.ply\\face_skin_hi.bmp");
    //M1.open_texture("/home/mepp/Desktop/face.texture.ply/face_skin_hi.bmp");
	
    M2.open_texture("C:\\Users\\noname\\Desktop\\face.obj\\face_skin_hi.jpg");
    //M2.open_texture("/home/mepp/Desktop/face.obj/face_skin_hi.jpg");

#if(0)
	OpenMesh::IO::Options opt;

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
    Build_polyhedron<Polyhedron1, HalfedgeDS1, Mesh1, MyMesh1> build_polyhedron1(M1, P1);
    P1.delegate(build_polyhedron1);
	std::cout << "--> (vertices: " << P1.size_of_vertices() << " - faces: " << P1.size_of_facets() << " - edges: " << (P1.size_of_halfedges() >> 1) << " - radius: " << P1.radius_ << " - n_scale: " << P1.normal_scale_ << ")\n\n";

	std::cout << "--> OpenMeshToCGALConverter\n";
	Polyhedron2 P2;
    Build_polyhedron<Polyhedron2, HalfedgeDS2, Mesh2, MyMesh2> build_polyhedron2(M2, P2);
    P2.delegate(build_polyhedron2);
	std::cout << "--> (vertices: " << P2.size_of_vertices() << " - faces: " << P2.size_of_facets() << " - edges: " << (P2.size_of_halfedges() >> 1) << " - radius: " << P2.radius_ << " - n_scale: " << P2.normal_scale_ << ")\n\n";

	// ---

	std::cout << "--> CGAL triangulation\n";
	P2.triangulate();
	std::cout << "--> (vertices: " << P2.size_of_vertices() << " - faces: " << P2.size_of_facets() << " - edges: " << (P2.size_of_halfedges() >> 1) << ")\n\n";

	// ---

	std::cout << "--> CGALToOpenMeshConverter\n";
	Mesh1 M1b;
	Build_mesh<Mesh1, MyMesh1, Polyhedron1> build_mesh1b(P1);
	build_mesh1b.build(M1b);
	std::cout << "--> (vertices: " << M1b.mesh().n_vertices() << " - faces: " << M1b.mesh().n_faces() << " - edges: " << M1b.mesh().n_edges() << " - radius: " << M1b.radius_ << " - n_scale: " << M1b.normal_scale_ << ")\n\n";

	std::cout << "--> CGALToOpenMeshConverter\n";
	Mesh2 M2b;
	Build_mesh<Mesh2, MyMesh2, Polyhedron2> build_mesh2b(P2);
	build_mesh2b.build(M2b);
	std::cout << "--> (vertices: " << M2b.mesh().n_vertices() << " - faces: " << M2b.mesh().n_faces() << " - edges: " << M2b.mesh().n_edges() << " - radius: " << M2b.radius_ << " - n_scale: " << M2b.normal_scale_ << ")\n\n";

	//return app.exec();

	// ---

	OpenMesh::IO::Options optw1;
	OpenMesh::IO::Options optw2;

	std::cout << "--> Writing mesh\n";
	if (M1b.mesh().has_vertex_colors()) { optw1 += OpenMesh::IO::Options::VertexColor; std::cout << "Mesh provides vertex colors\n"; }
	if (OpenMesh::IO::write_mesh(M1b.mesh(), "mesh_c.off", optw1))
		std::cout << "--> Writing done\n\n";
	else
		std::cout << "--> Writing error\n\n";

	std::cout << "--> Writing mesh\n";
	if (M2b.mesh().has_vertex_texcoords2D()) { optw2 += OpenMesh::IO::Options::VertexTexCoord; std::cout << "Mesh provides texture coordinates\n"; }
	if (OpenMesh::IO::write_mesh(M2b.mesh(), "mesh_t.off", optw2))
		std::cout << "--> Writing done\n\n";
	else
		std::cout << "--> Writing error\n\n";

    return app.exec();
}