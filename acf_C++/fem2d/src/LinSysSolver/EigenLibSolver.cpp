//
// EigenLibSolver.cpp
//
// acf_C++/fem2d
// Created by Wei Chen on 7/19/21
//

#include "EigenLibSolver.hpp"

namespace fem2d{

void EigenLibSolver::set_pattern(const std::vector<std::set<int>> &vNeighbor){
    Base::set_pattern(vNeighbor);

    coefMtr.resize(Base::numRows, Base::numRows);
    coefMtr.reserve(Base::nnz);
    
    memcpy(coefMtr.innerIndexPtr(), Base::ja.data(), Base::ja.size()*sizeof(Base::ja[0]));
    memcpy(coefMtr.outerIndexPtr(), Base::ia.data(), Base::ia.size()*sizeof(Base::ia[0]));
}

void EigenLibSolver::analyze_pattern(void){
    simplicialLDLT.analyzePattern(coefMtr);
    assert(simplicialLDLT.info() == Eigen::Success);
}

void EigenLibSolver::factorize(void){
    simplicialLDLT.factorize(coefMtr);
    assert(simplicialLDLT.info() == Eigen::Success);
}

void EigenLibSolver::solve(const Eigen::VectorXd &rhs, Eigen::VectorXd &result){
    result = simplicialLDLT.solve(rhs);
    assert(simplicialLDLT.info() == Eigen::Success);
}

void EigenLibSolver::setCoeff(int rowI, int colI, double val){
    assert(rowI < Base::numRows);
    const auto finder = Base::IJ2aI[rowI].find(colI);
    assert(finder != Base::IJ2aI[rowI].end());
    Base::a[finder->second] = val;
    coefMtr.valuePtr()[finder->second] = val;
}

void EigenLibSolver::addCoeff(int rowI, int colI, double val){
    assert(rowI < Base::numRows);
    const auto finder = Base::IJ2aI[rowI].find(colI);
    assert(finder != Base::IJ2aI[rowI].end());
    Base::a[finder->second] += val;
    coefMtr.valuePtr()[finder->second] += val;
}

void EigenLibSolver::setZero(void){
    Base::setZero();
    memcpy(coefMtr.valuePtr(), Base::a.data(), Base::a.size()*sizeof(Base::a[0]));
}

void EigenLibSolver::setUnit_row(int rowI){
    assert(rowI < Base::numRows);
    for(const auto &colIter : Base::IJ2aI[rowI]){
        double tmp = (colIter.first == rowI);
        Base::a[colIter.second] = tmp;
        coefMtr.valuePtr()[colIter.second] = tmp;
    }
}

void EigenLibSolver::setZero_col(int colI, const std::set<int> &rowVIs){
    assert(colI < Base::numRows);
    for(const auto &rowI : rowVIs){
        assert(rowI < Base::numRows);
        const auto finder = Base::IJ2aI[rowI].find(colI);
        if(finder != Base::IJ2aI[rowI].end()){
            Base::a[finder->second] = 0.0;
            coefMtr.valuePtr()[finder->second] = 0.0;
        }
    }
}

void EigenLibSolver::print_for_debug(void) const{
    std::cout << "**********" << std::endl;
    std::cout << coefMtr << std::endl;
    std::cout << "**********" << std::endl;
}

} // namespace fem2d