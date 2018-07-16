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

template<int surfdim, int spacedim, class VECTOR>
class Surface
{
public:
    Surface (Mapping_Q<surfdim,spacedim> &mapping,
                               Surface<surfdim,spacedim> &surface,
                               Triangulation<surfdim,spacedim> &triangulation,
                               unsigned int fe_degree,
                               double left_spectrum_bound,
                               double right_spectrum_bound,
                               double error_tolerance);
    ~LaplaceBeltramiAdaptivity();
    
    void run();
    
    
    
private:
    void setup_system();
    void assemble_system();
    void solve();
    void output_results();
    void pde_adapt();
    void geometry_adapt();
    
    
    
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
