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
    std::cout << "Choose Manifold:";
    std::cin >> surface_case;
    
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
    
    //Create Laplace Beltrami Function
    //LB function should do all the writing to folders, it will need surface name, fe_degree, and spectrum bounds in the folder name
    LaplaceBeltramiAdaptivity<surface_dimension, space_dimension, Vector<double>> laplace_beltrami_adapt(triangulation, fe_degree, left_spectrum_bound, right_spectrum_bound, error_tolerance);
    
    laplace_beltrami_adapt.run();
    
    //
}
