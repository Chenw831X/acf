//
// Mesh.cpp
// acf_C++/fem2d_heat
//
// Created by Wei Chen on 7/23/21
//

#include "Mesh.hpp"
#include "utils.hpp"

namespace fem2dHeat{

Mesh::Mesh(void){
    if(!Utils::readNodes("../input/coordinates.dat", V)){
        std::cout << "[fem2dHeat] --error-- fail in reading nodes" << std::endl;
        exit(-1);
    }
    if(!Utils::readTriangles("../input/elements3.dat", F)){
        std::cout << "[fem2dHeat] --error-- fail in reading triangles" << std::endl;
        exit(-1);
    }
    if(!Utils::readDBC("../input/dirichlet.dat", DBC)){
        std::cout << "[fem2dHeat] --error-- fail in reading dirichlet boundary conditions" << std::endl;
        exit(-1);
    }
    if(!Utils::readNBC("../input/neumann.dat", NBC)){
        std::cout << "[fem2dHeat] --error-- fail in reading nuemman boundary conditions" << std::endl;
        exit(-1);
    }
    
    std::cout << "[fem2dHeat] --info-- read mesh" << std::endl;
    // print_for_debug();

    computeFeatures();
    std::cout << "compute features" << std::endl;
}

Mesh::~Mesh(void){

}

void Mesh::print_for_debug(void){
    std::cout << "V=" << std::endl;
    std::cout << V << std::endl;
    std::cout << "F=" << std::endl;
    std::cout << F << std::endl;
    std::cout << "DBC=" << std::endl;
    std::cout << DBC << std::endl;
    std::cout << "NBC=" << std::endl;
    std::cout << NBC << std::endl;
}

void Mesh::computeFeatures(void){ // compute vNeighbor, FixedV, FreeV
    
}

} // namespace