//
// main.cpp
// acf_C++/fem2d
//
// Created by Wei Chen on 7/13/21
//

#include "Mesh.hpp"

#include <iostream>
#include <Eigen/Eigen>

int main(void){
    fem2d::Mesh mesh = fem2d::Mesh();

    mesh.assembly_stiffnessMat();
    std::cout << "Logging: [fem2d] assembly stiffness matrix" << std::endl;
    // mesh.linSysSolver->print_for_debug();

    mesh.assembly_rhs();
    std::cout << "Logging: [fem2d] assembly right hand side" << std::endl;

    mesh.adjustByDirichlet();
    std::cout << "Logging: [fem2d] adjust A and b based on dirichlet boundary conditions" << std::endl;

    mesh.solve();

    return 0;
}