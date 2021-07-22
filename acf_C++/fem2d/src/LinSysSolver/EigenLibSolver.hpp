//
// EigenLibSolver.hpp
//
// acf_C++/fem2d
// Created by Wei Chen on 7/19/21
//

#ifndef EigenLibSolver_hpp
#define EigenLibSolver_hpp

#include "LinSysSolver.hpp"

#include <Eigen/Eigen>
#include <vector>
#include <set>

namespace fem2d{

class EigenLibSolver : public LinSysSolver{
    typedef LinSysSolver Base;

protected:
    Eigen::SparseMatrix<double> coefMtr;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> simplicialLDLT;

public:
    virtual void set_pattern(const std::vector<std::set<int>> &vNeighbor);

    virtual void analyze_pattern(void);
    virtual void factorize(void);
    virtual void solve(const Eigen::VectorXd &rhs, Eigen::VectorXd &result);

    virtual void setCoeff(int rowI, int colI, double val);
    virtual void addCoeff(int rowI, int colI, double val);
    virtual void setZero(void);
    virtual void setUnit_row(int rowI);
    virtual void setZero_col(int colI, const std::set<int> &rowVIs);

    virtual void print_for_debug(void) const;
};

} // namespace

#endif // EigenLibSolver_hpp