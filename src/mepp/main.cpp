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
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>	MyMesh;

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


template <typename M>
bool 
open_mesh(const char* _filename, OpenMesh::IO::Options _opt)
{
	M mesh_;
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
	open_mesh<MyMesh>("C:\\Users\\noname\\Desktop\\face.texture.ply\\face.ply", opt);
	open_mesh<MyMesh2>("C:\\Users\\noname\\Desktop\\face.obj\\face.obj", opt);

	open_texture<MyMesh>("C:\\Users\\noname\\Desktop\\face.texture.ply\\face_skin_hi.bmp");
	open_texture<MyMesh2>("C:\\Users\\noname\\Desktop\\face.obj\\face_skin_hi.jpg");

    return app.exec();
}