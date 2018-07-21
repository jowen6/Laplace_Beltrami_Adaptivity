//
//  Surface.hpp
//  
//
//  Created by Justin Owen on 7/13/18.
//

#ifndef Surface_hpp
#define Surface_hpp

#include <stdio.h>


using namespace dealii;
//Set_Surface
//  attach manifold to triangulation
//  attach lift to triangulation
//  
template <int spacedim>
class LiftBase :  public Function<spacedim>
{
public:
    Lift () : Function<spacedim>(spacedim) {};
    
    virtual void vector_value (const Point<spacedim> &p,
                               Vector<double>   &values) const;
    
    
    virtual void vector_gradient(const Point<spacedim> &p, std::vector<Tensor<1,spacedim,double>> &gradients) const;
    
};


template<int surfdim, int spacedim, class VECTOR>
class Surface
{
public:
    Surface (Mapping_Q<surfdim,spacedim> &mapping,
                               Surface<surfdim,spacedim> &surface,
                               Triangulation<surfdim,spacedim> &triangulation,
                               unsigned int mapping_degree,

    ~Surface();
    
    
    
    
private:
    void lift();
    
    
    
    Triangulation<surfdim, spacedim>    &triangulation;
    FE_Q<surfdim, spacedim>             fe_basis;
    DoFHandler<surfdim, spacedim>       dof_handler;
    Surface<surfdim, spacedim>          &surface;
    
    
    
    ConstraintMatrix    stiffness_constraints;
    ConstraintMatrix    mass_constraints;
    SparsityPattern     sparsity_pattern;
    
    
    
    const std::string output_path;
    
};

#endif /* Surface_hpp */
