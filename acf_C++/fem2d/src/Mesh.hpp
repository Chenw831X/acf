//
// Mesh.hpp
// acf_C++/fem2d
//
// Created by Wei Chen on 7/9/21
//

#ifndef Mesh_hpp
#define Mesh_hpp

#include "EigenLibSolver.hpp"
#include <Eigen/Eigen>

#include <vector>
#include <set>

namespace fem2d{

class Mesh{
public: // owned data
    Eigen::MatrixXd V; // node coordinates
    Eigen::MatrixXi F3; // node indices of triangle elements
    Eigen::MatrixXi F4; // node indices of parallelograms elements
    Eigen::MatrixXi DBC; // Dirichlet boundary conditions
    Eigen::MatrixXi NBC; // Neumann boundary conditions

public: // owned features
    LinSysSolver *linSysSolver; // linear system solver
    std::vector<std::set<int>> vNeighbor; // records all vertice' indices adjacent to each vertice
    std::set<int> DBCV;
    std::set<int> FreeV;
    Eigen::VectorXd b;
    Eigen::VectorXd u;

public: // constructor
    Mesh(void);
    ~Mesh(void);

public:
    void print_for_debug(void) const;
    bool writeResult(void) const;
    void computeFeatures(void);

public: // API
    void assembly_stiffnessMat(void);
    void assembly_rhs(void);
    void adjustByDirichlet(void);
    void solve(void);

public: // interface: modify these 5 functions to define the problem
    void stima3(int triI, Eigen::Matrix3d &local);
    void stima4(int paraI, Eigen::Matrix4d &local);
    double f(const Eigen::Matrix<double, 1, 2> &Vmed);
    double g(const Eigen::Matrix<double, 1, 2> &Vmed);
    void u_d(void);
};

} // namespace fem2d

#endif // Mesh_hpp