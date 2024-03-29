/*!
 * \file mesh_to_polyhedron.h
 * \brief Class: Build_polyhedron and Build_mesh.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2013
 */

#ifndef _MESH_TO_POLYHEDRON_
#define _MESH_TO_POLYHEDRON_

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

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <OpenMesh/Core/IO/importer/ImporterT.hh>

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

// A modifier creating a polyhedron with the incremental builder.
template <class P, class HDS, class Mesh, class MyMesh>
class Build_polyhedron : public CGAL::Modifier_base<HDS>
{
public:
    Build_polyhedron(Mesh& _mesh, P& _p) : mesh_(_mesh.mesh()), p_(_p) { p_.tex_id(_mesh.tex_id()); }

    void operator()(HDS& hds)
    {
		// Postcondition: 'hds' is a valid polyhedral surface.
		CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

		B.begin_surface(mesh_.n_vertices(), mesh_.n_faces(), mesh_.n_halfedges());

        typedef typename HDS::Vertex                    Vertex;
        typedef typename Vertex::Point                  Point;

        typedef typename HDS::Vertex_handle             Vertex_handle;
        typedef typename HDS::Face_handle               Face_handle;

        // vertex
        typename MyMesh::VertexIter vit;
        for(vit=mesh_.vertices_begin(); vit!=mesh_.vertices_end(); ++vit)
        {
                Vertex_handle vertex = B.add_vertex(Point(mesh_.point(vit)[0], mesh_.point(vit)[1], mesh_.point(vit)[2]));

                if (mesh_.has_vertex_colors())
                {
                        OpenMesh::Vec3uc uc = mesh_.color(vit);
                        vertex->color(/*(float)*/uc[0]/*/255.0*/, /*(float)*/uc[1]/*/255.0*/, /*(float)*/uc[2]/*/255.0*/); // TODO alpha (colorA)
                }

				if (mesh_.has_vertex_texcoords2D())
                {
                        OpenMesh::Vec2f tc = mesh_.texcoord2D(vit);
                        vertex->texcoord2D(tc[0], tc[1]);
                }
        }

        // facet
        typename MyMesh::ConstFaceIter          fIt(mesh_.faces_begin()), fEnd(mesh_.faces_end());
        typename MyMesh::ConstFaceVertexIter    fvIt;
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
                        face->color(/*(float)*/uc[0]/*/255.0*/, /*(float)*/uc[1]/*/255.0*/, /*(float)*/uc[2]/*/255.0*/); // TODO alpha (colorA)
                }
        }

		B.end_surface();

        if (B.check_unconnected_vertices())
                B.remove_unconnected_vertices();

        // ---

		// normals
        p_.compute_normals();

		if (p_.size_of_vertices()>0)
		{
			// colors
			if (mesh_.has_vertex_colors())
					p_.has_vertex_colors(true);
			if (mesh_.has_face_colors())
					p_.has_face_colors(true);
			if (p_.has_vertex_colors())
					std::cout << "Polyhedron provides vertex colors\n";
			if (p_.has_face_colors())
					std::cout << "Polyhedron provides face colors\n";

			// texture
			if (mesh_.has_vertex_texcoords2D())
					p_.has_vertex_texcoords2D(true);
			if (p_.has_vertex_texcoords2D())
					std::cout << "Polyhedron provides texture coordinates\n";

			// bounding box
			p_.compute_bounding_box();

			// set center and radius ( set_scene_pos( , ) )
			p_.center( (p_.bbMin() + p_.bbMax())*0.5 );
			p_.radius( std::sqrt( ::CGAL::to_double( (p_.bbMin() - p_.bbMax()).squared_length() ) )*0.5 );
                
			// set normal scale for normal display
			typename P::Vector diff = p_.bbMax() - p_.bbMin();
			p_.normal_scale( ::CGAL::to_double( ::CGAL::min( ::CGAL::min(diff[0], diff [1]), diff[2] ) )*0.05f );

			// TODO set base point for displaying face normals
		}
	}

private:
	MyMesh& mesh_;
	P& p_;
};

template <class Mesh, class MyMesh, class P>
class Build_mesh
{
public:
	Build_mesh(P& _p) : p_(_p) {}

	void build(Mesh& _mesh)
    {
            MyMesh& mesh_ = _mesh.mesh();

			mesh_.request_face_normals();
			mesh_.request_face_colors();
			mesh_.request_vertex_normals();
			mesh_.request_vertex_colors();
			mesh_.request_vertex_texcoords2D();
                
            typedef typename P::Vertex_iterator						VI;
            typedef typename P::Facet_iterator						FI;
            typedef typename P::Halfedge_around_facet_circulator	HFC;

            OpenMesh::IO::ImporterT<MyMesh> importerOpenMesh(mesh_);
            importerOpenMesh.reserve(p_.size_of_vertices(), 3*p_.size_of_vertices(), p_.size_of_facets());

            OpenMesh::Vec3f v;
            OpenMesh::VertexHandle vh;

            OpenMesh::IO::BaseImporter::VHandles vhandles;
            OpenMesh::FaceHandle fh;

            // output vertices
            for (VI vi = p_.vertices_begin(); vi != p_.vertices_end(); ++vi)
            {
                    v[0] = ::CGAL::to_double(vi->point().x());
                    v[1] = ::CGAL::to_double(vi->point().y());
                    v[2] = ::CGAL::to_double(vi->point().z());

                    vh = importerOpenMesh.add_vertex(v);

                    if (p_.has_vertex_colors())
                            mesh_.set_color(vh, OpenMesh::Vec3uc(vi->color(0)/**255*/, vi->color(1)/**255*/, vi->color(2)/**255*/)); // TODO alpha (colorA)

					if (p_.has_vertex_texcoords2D())
                            mesh_.set_texcoord2D(vh, OpenMesh::Vec2f(vi->texcoord2D(0), vi->texcoord2D(1)));
            }

            // constructs an inverse index
            typedef CGAL::Inverse_index<VI> Index;
            Index index(p_.vertices_begin(), p_.vertices_end());

            // output facets
            for (FI fi = p_.facets_begin(); fi != p_.facets_end(); ++fi)
            {
                    vhandles.clear(); // OpenMesh

                    HFC hc = fi->facet_begin();
                    HFC hc_end = hc;
                    std::size_t n = circulator_size(hc);
                    CGAL_assertion(n >= 3);

                    // facet begin
                    do
                    {
                            vhandles.push_back(OpenMesh::VertexHandle(index[VI(hc->vertex())]));
                            ++hc;
                    } while (hc != hc_end);
                    // facet end

                    fh = importerOpenMesh.add_face(vhandles);

                    if (p_.has_face_colors())
                            mesh_.set_color(fh, OpenMesh::Vec3uc(fi->color(0)/**255*/, fi->color(1)/**255*/, fi->color(2)/**255*/)); // TODO alpha (colorA)
            }
                
            // ---

			// normals
			mesh_.update_face_normals();
            mesh_.update_vertex_normals();

			if (mesh_.n_vertices()>0)
			{
				// colors
				if (!p_.has_vertex_colors())
					mesh_.release_vertex_colors();
				if (!p_.has_face_colors())
					mesh_.release_face_colors();
				if (mesh_.has_vertex_colors())
					std::cout << "Mesh provides vertex colors\n";
				if (mesh_.has_face_colors())
					std::cout << "Mesh provides face colors\n";
			
				// texture
				if (!p_.has_vertex_texcoords2D())
					mesh_.release_vertex_texcoords2D();
				if (mesh_.has_vertex_texcoords2D())
					std::cout << "Mesh provides texture coordinates\n";
				_mesh.tex_id(p_.tex_id());

				// bounding box
				typename MyMesh::ConstVertexIter vIt(mesh_.vertices_begin());
				typename MyMesh::ConstVertexIter vEnd(mesh_.vertices_end());

				_mesh.bbMin( OpenMesh::vector_cast<OpenMesh::Vec3f>(mesh_.point(vIt)) ); _mesh.bbMax( _mesh.bbMin() );
                
				for (size_t count=0; vIt!=vEnd; ++vIt, ++count)
				{
					_mesh.bbMin().minimize( OpenMesh::vector_cast<OpenMesh::Vec3f>(mesh_.point(vIt)) );
					_mesh.bbMax().maximize( OpenMesh::vector_cast<OpenMesh::Vec3f>(mesh_.point(vIt)) );
				}
                    
				// set center and radius ( set_scene_pos( , ) )
				_mesh.center( (_mesh.bbMin() + _mesh.bbMax())*0.5 );
				_mesh.radius( (_mesh.bbMin() - _mesh.bbMax()).norm()*0.5 );
                
				// set normal scale for normal display
				_mesh.normal_scale( (_mesh.bbMax() - _mesh.bbMin()).min()*0.05f );

				// TODO set base point for displaying face normals
			}
	}

private:
	P& p_;
};

#endif // _MESH_TO_POLYHEDRON_
