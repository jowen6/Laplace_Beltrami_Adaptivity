//
//  Laplace_Beltrami_Adapt.hpp
//  
//
//  Created by Justin Owen on 7/13/18.
//

#ifndef Laplace_Beltrami_Adapt_hpp
#define Laplace_Beltrami_Adapt_hpp

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;


template<int surfdim, int spacedim, class VECTOR>
class LaplaceBeltramiAdaptivity
{
    public:
    LaplaceBeltramiAdaptivity (const std::string output_path,
                               Triangulation<surfdim,spacedim> &triangulation,
                               MappingQ<surfdim, spacedim> &mapping,
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
    MappingQ<surfdim, spacedim>         &mapping;
    
    
    
    ConstraintMatrix    stiffness_constraints;
    ConstraintMatrix    mass_constraints;
    SparsityPattern     sparsity_pattern;
    
    
    
    const std::string output_path;
    
};

#endif /* Laplace_Beltrami_Adapt_hpp */
