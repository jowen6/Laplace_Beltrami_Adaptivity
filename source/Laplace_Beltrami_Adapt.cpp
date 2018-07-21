//
//  Laplace_Beltrami_Adapt.cpp
//  
//
//  Created by Justin Owen on 7/13/18.
//

#include "Laplace_Beltrami_Adapt.hpp"
//Complete
template <int surfdim, int spacedim, class VECTOR>
LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::LaplaceBeltramiAdaptivity(const std::string output_path,
                                                                                Triangulation<surfdim, spacedim> &triangulation,
                                                                                MappingQ<surfdim, spacedim> &mapping
                                                                                unsigned int fe_degree,
                                                                                double left_spectrum_bound,
                                                                                double right_spectrum_bound,
                                                                                double error_tolerance):

    triangulation(triangulation),
    fe_basis(fe_degree),
    dof_handler(triangulation),
    lift_dof_handler(triangulation),
    mapping(mapping),
    output_path(output_path)




//Complete
template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::~LaplaceBeltramiAdaptivity(){
    
    lift_dof_handler.clear();
    dof_handler.clear();
    
}



//Incomplete
template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::setup_system(){

    dof_handler.distribute_dofs(fe_basis);
}



//Incomplete
template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::assemble_system(){
    
    stiff_matrix = 0;
    mass_matrix = 0;
    
    //Create quadrature rule
    unsigned int iquad = 2*(fe_basis.degree)+(mapping.get_degree) + (mapping.get_degree)*spacedim+1;
    QGauss<surfdim>     quadrature_formula(iquad);
    const unsigned int   n_q_points    = quadrature_formula.size();
    
    //Evaluate FEM basis functions at quadrature points
    FEValues<surfdim, spacedim>     fe_values(mapping,
                                              fe_basis,
                                              quadrature_formula,
                                              update_values |
                                              update_gradients |
                                              update_quadrature_points |
                                              update_JxW_values);
    
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;

    //Initialize local cell matrices
    FullMatrix<double>              cell_stiff_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double>              cell_mass_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<unsigneed int>      local_dof_indices(dofs_per_cell);
    
    
    typename DoFHandler<surfdim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    
    //Loop through elements
    for (; cell!=endc; ++cell){
        
        cell_stiff_matrix = 0;
        cell_mass_matrix = 0;
        
        fe_values.reinit(cell);
        
        for (unsigned int i=0; i<dofs_per_cell; ++i){
            for (unsigned int j=0; j<dofs_per_cell; ++j){
                for (unsigned int q_point=0; q_point<n_q_points; ++q_point){
                    
                    cell_stiff_matrix(i,j) += fe_values.shape_grad(i,q_point)
                                              * fe_values.shape_grad(j,q_point)
                                              * fe_values.JxW(q_point);
                    
                    cell_mass_matrix(i,j) += fe_values.shape_value(i,q_point)
                                             * fe_values.shape_value(j,q_point)
                                             * fe_values.JxW(q_point);

                }
            }
        }
    }
}



//Incomplete
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



//Incomplete
template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::pde_adapt() const{
    
    //Check PDE error
    setup_system();
    assemble_system();
    solve();
    double pde_estimate;
    //PDE Estimator
    
    while (pde_estimate > error_tolerance){
        

        GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                         estimated_error_per_cell,
                                                         0.3, 0.03);
    
        triangulation.execute_coarsening_and_refinement ();
        
        setup_system();
        assemble_system();
        solve();
        double pde_estimate;
        
    }
}


//Incomplete
template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::approximate_lift() const{
    
    //Take lift function and interpolate to create approximate_lift
    
    
}


//Incomplete
template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::geometry_adapt() const{
    
    //error_tolerance fraction
    double omega = 0.3;

    //Geometric estimator on each cell
    VectorTools::integrate_difference (mapping,
                                       approximate_lift,
                                       lift,
                                       estimator_geometric_error_per_cell,
                                       QGauss<2>(mapping_degree+1),
                                       VectorTools::W1infty_seminorm);
    
    //Geometric estimate for loop
    double geometric_estimate = estimator_geometric_error_per_cell.linfty_norm();
    
    //Estimate and refine until geometric error is sufficiently small
    while (omega*error_tolerance < geometric_estimate){
        
        //Input vector from estimator and greedy marks 30% of total error
        GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                         estimator_geometric_error_per_cell,
                                                         0.3, 0.03);
        
        //Refine triangulation
        triangulation.execute_coarsening_and_refinement ();
        
        //Geometric estimator on each cell
        VectorTools::integrate_difference (mapping,
                                           approximate_lift,
                                           lift,
                                           estimator_geometric_error_per_cell,
                                           QGauss<2>(mapping_degree+1),
                                           VectorTools::W1infty_seminorm);
        
        //Geometric estimate for loop
        geometric_estimate = estimator_geometric_error_per_cell.linfty_norm();
       
        
    }
}



//Incomplete
template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::output_results() const{
    std::string folder = output_path
                         +"_FEMDeg" + std::to_string(fe_degree)
                         +"_LSpec" + left_spectrum
                         +"_RSpec" + right_spectrum;
    
    MyUtilities::mkpath(folder.c_str(),0777);
    
    std::string eigenvalues_filename = folder + "/eigenvalues.txt";
    std::remove(eigenvalues_filename.c_str());
    
    std::string eigenvectors_filename = folder + "/eigenvectors.vtk";
    std::remove(eigenvectors_filename.c_str());
    
    std::string pde_estimator_filename = folder + "/pde_error_estimator.txt";
    std::remove(pde_estimator_filename.c_str());
    
    std::string geo_estimator_filename = folder + "/geo_error_estimator.txt";
    std::remove(geo_estimator_filename.c_str());
    
    
    //Output Eigenvalues
    std::ofstream of_eigenvalues(eigenvalues_filename,std::ios::app);
    of_eigenvalues.precision(10);
    of_eigenvalues<<std::scientific;
    
    for (int i=0;i<num_eigenvalues;++i){
        of_eigenvalues<<eigenvalues[i]<<std::endl;
    }
    
    
    //Initialize data output for eigenvectors
    DataOut<surfdim,DoFHandler<surfdim,spacedim> > data_out;
    data_out.attach_dof_handler (dof_handler);
    
    //write eigenvectors to data_out
    for (unsigned int i=0; i<eigenvectors.size(); ++i){
        localized_solutions[i].reinit(dof_handler.n_dofs());
        localized_solutions[i]=eigenfunctions[i];
        
        data_out.add_data_vector (localized_solutions[i],
                                  std::string("eigenvector_") + Utilities::int_to_string(i),
                                  DataOut<surfdim, DoFHandler<surfdim,spacedim> >::type_dof_data);
    }
    
    data_out.build_patches (mapping, mapping.get_degree());
    
    //Output Eigenvectors
    std::ofstream output (eigenvectors_filename);
    data_out.write_vtk (output);


}



template <int surfdim, int spacedim, class VECTOR>
void LaplaceBeltramiAdaptivity<surfdim, spacedim, VECTOR>::run(){
    
    geometry_adapt();
    pde_adapt();
    output_results();
    
}
