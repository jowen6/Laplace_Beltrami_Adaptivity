//
//  Laplace_Beltrami_Adapt.cpp
//  
//
//  Created by Justin Owen on 7/13/18.
//

#include "Laplace_Beltrami_Adapt.hpp"

template <int surfdim, int spacedim, class VECTOR>
LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::LaplaceBeltramiAdaptivity(const std::string output_path,
                                                                                Triangulation<surfdim, spacedim> &triangulation,
                                                                                unsigned int fe_degree,
                                                                                double left_spectrum_bound,
                                                                                double right_spectrum_bound,
                                                                                double error_tolerance):

    triangulation(triangulation),
    fe_basis(fe_degree),
    dof_handler(triangulation),
    mapping(mapping),
    output_path(output_path)





template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::~LaplaceBeltramiAdaptivity(){
    
    dof_handler.clear();
}



template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::setup_system(){

    dof_handler.distribute_dofs(fe_basis);
}



template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::assemble_system(){
    
    
}



template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::solve() const{
    
    SolverControl solver_control (2*dh.n_dofs(), 1e-10,true,true);
    unsigned int fe_degree = fe.degree;
    unsigned int mapping_degree = mapping.get_degree();
    double mean;
    
    mpi_mgr.write("Finding the Eigenvalues");
    
    //SolverControl solver_control (dh.n_dofs(), 1e-10,true,true);
    SLEPcWrappers::SolverKrylovSchur eigensolver(solver_control,mpi_mgr.get_communicator());
    
    //eigensolver.set_which_eigenpairs (EPS_SMALLEST_REAL); //Find eigenvalues with smallest real components
    eigensolver.set_problem_type (EPS_GHEP);  //If not specified, eigenvalues with multiplicity >1 share a single eigenvector
    //Doesn't seem like SLEPc can handle solving Hermitian system if it assumes non-Hermitian
    
    //override some options
    //    EPSSetType(*eigensolver.get_EPS(),EPSKRYLOVSCHUR);
    
    
    eigensolver.solve (system_matrix, mass_matrix,
                       eigenvalues, eigenfunctions,
                       eigenfunctions.size());
    
    
    PETScWrappers::Vector  localized_solution;  // not useful anymore
    //    double mean;
    //    unsigned int fe_degree = fe.degree;
    //    unsigned int mapping_degree = mapping.get_degree();
    
    
    localized_solution.reinit(dh.n_dofs());
    
    for (unsigned int i=0; i<eigenfunctions.size(); ++i){
        localized_solution=eigenfunctions[i];
        // normalization  this does not work in parallel, use MPI_REDUCE
        mean = VectorTools::compute_mean_value(mapping, dh, QGauss<dim>(2*(fe_degree)+(mapping_degree)+(mapping_degree)*dim+1), localized_solution, 0);
        // substract off mean value
        for (unsigned int i = 0;i<dh.n_dofs();++i)
            localized_solution(i)-=mean;
        
        localized_solution.compress(VectorOperation::add);
        
        eigenfunctions[i]=localized_solution;
        eigenfunctions[i].compress(VectorOperation::insert);
        
    }
    
    
}



template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::pde_adapt() const{
    
    //geometry_adapt has already been run
    //while eta>epsilon
    //setup_system
    //assemble system
    //solve
    //PDE Estimator
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     estimated_error_per_cell,
                                                     0.3, 0.03);
    
    triangulation.execute_coarsening_and_refinement ();
    
}

template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::geometry_adapt() const{

    //Expecting a number of initial global refinements of triangulation
    //while zeta>omega*epsilon
    //GEO Estimator
    //Input vector from estimator and greedy marks 30% of total error
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     estimated_error_per_cell,
                                                     0.3, 0.03);
    //Refine triangulation
    triangulation.execute_coarsening_and_refinement ();
    
    
}


template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::output_results() const{
    

    
}



template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::run(){
    
    geometry_adapt();
    pde_adapt();
    output_results();
    
}
