//
// EigenLibSolver.hpp
// acf_C++/fem2d_heat
//
// Created by Wei Chen on 7/24/21
//

#ifndef EigenLibSolver_hpp
#define EigenLibSolver_hpp

#include "LinSysSolver.hpp"

namespace fem2dHeat{

class EigenLibSolver : public LinSysSolver{
    typedef LinSysSolver Base;

protected: // owned features
    Eigen::SparseMatrix<double> coefMtr;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> simplicialLDLT;

public:
    EigenLibSolver operator + (const EigenLibSolver &other);

public: // API
    virtual void set_pattern(const std::vector<std::set<int>> &vNeighbor);

    virtual void analyze_pattern(void);
    virtual void factorize(void);
    virtual void solve(const Eigen::VectorXd &rhs, Eigen::VectorXd &res);

    virtual void setCoeff(int rowI, int colI, double val);
    virtual void addCoeff(int rowI, int colI, double val);
    virtual void setZero(void);
    virtual void setUnit_row(int rowI);
    virtual void setZero_col(int colI, const std::set<int> &rowVIs);
};

} // namespace

#endif