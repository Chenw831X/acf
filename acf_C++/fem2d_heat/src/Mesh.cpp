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
    std::cout << "[fem2dHeat] --info-- compute features" << std::endl;

    linSysSolver = new EigenLibSolver();
    linSysSolver->set_pattern(vNeighbor);
    std::cout << "[fem2dHeat] --info-- linear solver: set pattern" << std::endl;
}

Mesh::~Mesh(void){
    delete linSysSolver;
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
    vNeighbor.resize(V.rows()); // compute vNeighbor
    for(int vI=0; vI<V.rows(); ++vI){
        vNeighbor[vI].insert(vI);
    }
    for(int triI=0; triI<F.rows(); ++triI){
        const Eigen::Matrix<int, 1, 3> &triVInd = F.row(triI);
        for(int i=0; i<3; ++i){
            int vI = triVInd(i);
            for(int j=i+1; j<3; ++j){
                int vJ = triVInd(j);
                vNeighbor[vI].insert(vJ);
                vNeighbor[vJ].insert(vI);
            }
        }
    }

    for(int DBCI=0; DBCI<DBC.rows(); ++DBCI){ // compute FixedV
        const Eigen::Matrix<int, 1, 2> &DBCVInd = DBC.row(DBCI);
        FixedV.insert(DBCVInd(0));
        FixedV.insert(DBCVInd(1));
    }

    for(int vI=0; vI<V.rows(); ++vI){ // compute FreeV
        if(FixedV.find(vI)==FixedV.end()){
            FreeV.insert(vI);
        }
    }
}

void Mesh::assembly_K(void){
    
}

    /*
    void assembly_b(void);
    void adjustByDirichlet(void);
    void solve(void);
    */
} // namespace