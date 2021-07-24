//
// Mesh.hpp
// acf_C++/fem2d_heat
//
// Created by Wei Chen on 7/23/21
//

#ifndef Mesh_hpp
#define Mesh_hpp

#include "EigenLibSolver.hpp"

#include <Eigen/Eigen>

#include <iostream>
#include <vector>
#include <set>
#include <map>

namespace fem2dHeat{

class Mesh{
public: // owned data
    Eigen::MatrixXd V; // coordinates of mesh
    Eigen::MatrixXi F; // triangle elements of mesh
    Eigen::MatrixXi DBC; // Dirichlet boundary conditions, each row contains two indices of an edge
    Eigen::MatrixXi NBC; // Neumman boundary conditions, each row contains two indices of an edge

public: // owned features
    std::vector<std::set<int>> vNeighbor;
    std::set<int> FixedV;
    std::set<int> FreeV;
    LinSysSolver *linSysSolver;
    
public: // constructor
    Mesh(void);
    ~Mesh(void);

public:
    void computeFeatures(void);

public: // API
    void print_for_debug(void);
    void assembly_K(void);
    void assembly_b(void);
    void adjustByDirichlet(void);
    void solve(void);

public: // interfaces

};

} // namespace

#endif // Mesh_hpp