//
// main.cpp
// acf_C++/fem2d_heat
//
// Created by Wei Chen on 7/23/21
//

#include "Mesh.hpp"

#include <Eigen/Eigen>

#include <iostream>

int main(){
    fem2dHeat::Mesh mesh = fem2dHeat::Mesh();
    std::cout << "[fem2dHeat] --info-- mesh constructed" << std::endl;

    return 0;
}