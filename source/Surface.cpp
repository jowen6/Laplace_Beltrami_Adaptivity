//
//  Surface.cpp
//  
//
//  Created by Justin Owen on 7/13/18.
//

#include "Surface.hpp"
template <int spacedim>
class LiftBase :  public Function<spacedim>
{
public:
    Lift () : Function<spacedim>(spacedim) {};
    
    virtual void vector_value (const Point<spacedim> &p,
                               Vector<double>   &values) const;
    
    
    virtual void vector_gradient(const Point<spacedim> &p, std::vector<Tensor<1,spacedim,double>> &gradients) const;
    
};
/*
 template <int surfdim, spacedim>
 class SurfaceBase {
 virtual void Set_geometry(triangulation);
 Liftbase<spacedim> Lift();
 
 };
 
 template <int surfdim, int spacedim>
 class Sphere : public SurfaceBase <spacedim> {
 virtual void Set_geometry(triangulation);
 Liftbase<spacedim> Lift();
 
 };
*/
 

template<int surfdim, int spacedim>
void Surface(int CASE, Triangulation<surfdim,spacedim> &triangulation)
{
    switch (CASE){

        case 0: // SPHERE
        {
            /*
            static SphericalManifold<2,3> surface_description;
            GridGenerator::hyper_sphere(triangulation);
            triangulation.set_all_manifold_ids(0);
            triangulation.set_manifold (0, surface_description);
            */
            
            static SphericalManifold<2,3> sphere_surface_description;
            static FlatManifold<2,3> polyhedral_surface_description;
            GridGenerator::hyper_sphere (triangulation);
            triangulation.set_all_manifold_ids(0);
            triangulation.set_manifold(0, sphere_surface_description);
            
            //triangulation.refine_global(5);
            
            triangulation.refine_global(2);
            triangulation.set_manifold(0, polyhedral_surface_description);
            
            break;
        }

        case 1: // Pi/3 sphere
        {
            static SphericalManifold<2,3> surface_description;
            {
                Triangulation<3> volume_mesh;
                GridGenerator::half_hyper_ball(volume_mesh);
                std::set<types::boundary_id> boundary_ids;
                boundary_ids.insert (0);
                GridGenerator::extract_boundary_mesh (volume_mesh, triangulation,
                                                      boundary_ids);
            }
            triangulation.set_all_manifold_ids(0);
            triangulation.set_manifold (0, surface_description);

            break;
        }
            
        case 2: // triangle
        {
            triangle_builder.buildTriangulation(&triangulation,
                                                Point<3>(0.,0.,0.),
                                                Point<3>(4.0,0.,0.),
                                                Point<3>(0.,3.0,0.));
            
            break;
        }
            
        case 3: // torus
        {
            GridGenerator::torus(triangulation,Rout,Rin);
            triangulation.set_boundary(0,torus_boundary_description);
            
            break;
        }
            
        case 4: // glowinski
            // requires c_p1, c_p2, c_p2 and c_coeff
            // the surface is c_coeff x y z (1-x/c_p1-y/c_p2-z/c_p3 = 1
        {
            
            if (deal_II_dimension != 3){
                mpi_mgr.write("Wrong dimension");
                exit(0);
            }
            
            c_coeff= 50.0;
            c_p1 = 2.0;
            c_p2 = 3.0;
            c_p3 = 5.0;
            
            tetrahedra_builder.buildTriangulation(&triangulation,
                                                  Point<3>(0.0,0.0,0.0),
                                                  Point<3>(c_p1,0.,0.),
                                                  Point<3>(0.,c_p2,0.),
                                                  Point<3>(0.,0.0,c_p3));
            
            static ImplicitManifold glowinski_manifold(c_coeff,c_p1,c_p2,c_p3);
            // project current vertices... yeah need to bind the class function
            std::function<Point<3>(const Point<3> &)> project_from_glowinski;
            
            project_from_glowinski = std::bind(&ImplicitManifold::project_to_manifold_divide,
                                               &glowinski_manifold,
                                               std::placeholders::_1);

            GridTools::transform (project_from_glowinski, triangulation);

            triangulation.set_all_manifold_ids(0);
            triangulation.set_manifold (0, glowinski_manifold);
            
            break;
        }
            
        case 5: // glowinski sym
            // requires c_p1, c_p2, c_p2 and c_coeff
            // the surface is c_coeff x y z (1-x/c_p1-y/c_p2-z/c_p3 = 1
        {
            
            if (deal_II_dimension != 3){
                mpi_mgr.write("Wrong dimension");
                exit(0);
            }
            
            c_coeff= 1024.0;
            c_p1 = 1.0;
            c_p2 = 1.0;
            c_p3 = 1.0;
            
            tetrahedra_builder.buildTriangulation(&triangulation,
                                                  Point<3>(0.,0.,0.),
                                                  Point<3>(c_p1,0.,0.),
                                                  Point<3>(0.,c_p2,0.),
                                                  Point<3>(0.,0.0,c_p3));
            
            static ImplicitManifold glowinski_manifold(c_coeff,c_p1,c_p2,c_p3);
            triangulation.set_all_manifold_ids(0);
            
            // project current vertices... yeah need to bind the class function
            std::function<Point<3>(const Point<3> &)> project_from_glowinski;
            
            project_from_glowinski = std::bind(&ImplicitManifold::project_to_manifold_divide,
                                               &glowinski_manifold,
                                               std::placeholders::_1);

            GridTools::transform (project_from_glowinski, triangulation);
            triangulation.set_manifold (0, glowinski_manifold);

            break;
        }

        case 6: // HALF SPHERE
            /*
             {
             GridGenerator::half_hyper_ball(volume_mesh);
             
             volume_mesh.set_boundary (1,hyper_ball_boundary_description);
             volume_mesh.set_boundary (0, hyper_ball_boundary_description);
             tria.set_boundary (0, surface_description);
             boundary_ids.insert(0);
             GridGenerator::extract_boundary_mesh (volume_mesh, tria,boundary_ids);
             }
             break;
             */
            
        case 7: // Rectangle
        {
            quadrilateral_builder.buildTriangulation(&triangulation,Point<3>(0,0,0),
                                                     Point<3>(1,0,0),
                                                     Point<3>(0,1.000,0),
                                                     Point<3>(1,1.0,0));

            triangulation.set_boundary(0);
            
            break;
            
        }
            
        case 8:
        {
            static SphericalManifold<2,3> surface_description;
            GridGenerator::hyper_sphere(triangulation);
            triangulation.set_all_manifold_ids(0);
            triangulation.set_manifold (0, surface_description);
            
            break;
            
        }
            
        case 9:
        {
            static SphericalManifold<2,3> surface_description;
            GridGenerator::hyper_sphere(triangulation);
            triangulation.set_all_manifold_ids(0);
            triangulation.set_manifold (0, surface_description);
            
            break;
        }
            
        case 10: // Dziuk Modified
        {
            GridGenerator::hyper_sphere(triangulation,Point<3>(0,0,0),sqrt(3));
            triangulation.set_boundary(0);

            static ImplicitManifoldDziuk dziuk_manifold;
            
            
            // project current vertices... yeah need to bind the class function
            std::function<Point<3>(const Point<3> &)> project_from_dziuk;
            project_from_dziuk = std::bind(&ImplicitManifoldDziuk::project_to_manifold_divide,
                                           &dziuk_manifold,
                                           std::placeholders::_1);

            GridTools::transform (project_from_dziuk, triangulation);
            
            triangulation.set_all_manifold_ids(0);
            triangulation.set_manifold (0, dziuk_manifold);
            
            break;
        }
            
        default:
            std::cout<<"CASE NOT IMPLEMENTED"<<endl;
            exit(1);
            
        }
}
