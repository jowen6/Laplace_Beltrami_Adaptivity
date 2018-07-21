//
//  main.cpp
//  Laplace_Beltrami_Adaptivity
//
//  Created by Justin Owen on 7/13/18.
//

#include "main.hpp"
#include "Laplace_Beltrami_Adapt.hpp"unsigned int fe_degree;

unsigned int init_refinement;
unsigned int surface_dimension = 2;
unsigned int space_dimension = 3;
double error_tolerance;
std::string left_spectrum_bound;
std::string right_spectrum_bound;
std::string surface_name;


int main ( int argc, char **argv ){
    
    //Ask user to specify parameters
    std::cout << "Choose Manifold:";
    std::cin >> surface_case;
    
    std::cout << "Enter Error Tolerance:";
    std:: cin >> error_tolerance;
    
    std::cout << "Enter Finite Element Degree:";
    std:: cin >> fe_degree;
    
    std::cout << "Initial Refinement Level:";
    std::cin >> init_refinement;
    
    std::cout << "Left Spectrum Bound:";
    std::cin >> left_spectrum_bound;
    
    std::cout << "Right Spectrum Bound";
    std::cin >> right_spectrum_bound;

    //Initialize Triangulation
    Triangulation<surface_dimension, space_dimension> triangulation;
    
    //Create Surface
    Surface<surface_dimension, space_dimension>(surface_case, &triangulation);
    
    //Refine Surface
    triangulation.refine_global(init_refinement);
    
    //Beginning of output path
    std::string output_folder = "output/Surface"+std::to_string(CASE);
    
    //Create Laplace Beltrami Class
    LaplaceBeltramiAdaptivity<surface_dimension, space_dimension, Vector<double>> laplace_beltrami_adapt(output_folder,
                               triangulation,
                               mapping,
                               fe_degree,
                               atoi(left_spectrum_bound.c_str()),
                               atoi(right_spectrum_bound.c_str()),
                               atoi(error_tolerance.c_str()));
    
    //Run Laplace Beltrami Code
    laplace_beltrami_adapt.run();
    
    //
}
