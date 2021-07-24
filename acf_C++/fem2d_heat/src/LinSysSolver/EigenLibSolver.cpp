//
// EigenLibSolver.cpp
// acf_C++/fem2d_heat
//
// Created by Wei Chen on 7/24/21
//

#include "EigenLibSolver.hpp"

namespace fem2dHeat{

EigenLibSolver EigenLibSolver::operator+(const EigenLibSolver &other){ // two operators must have same pattern
    
}

void EigenLibSolver::set_pattern(const std::vector<std::set<int>> &vNeighbor){
    Base::set_pattern(vNeighbor);

    coefMtr.resize(Base::numRows, Base::numRows);
    coefMtr.reserve(Base::nnz);

    memcpy(coefMtr.innerIndexPtr(), Base::ja.data(), Base::ja.size()*sizeof(Base::ja(0)));
    memcpy(coefMtr.outerIndexPtr(), Base::ia.data(), Base::ia.size()*sizeof(Base::ia(0)));
}
    
void EigenLibSolver::analyze_pattern(void){
    simplicialLDLT.analyzePattern(coefMtr);
    assert(simplicialLDLT.info() == Eigen::Success);
}

void EigenLibSolver::factorize(void){
    simplicialLDLT.factorize(coefMtr);
    assert(simplicialLDLT.info() == Eigen::Success);
}

void EigenLibSolver::solve(const Eigen::VectorXd &rhs, Eigen::VectorXd &res){
    res = simplicialLDLT.solve(rhs);
    assert(simplicialLDLT.info() == Eigen::Success);
}

void EigenLibSolver::setCoeff(int rowI, int colI, double val){
    assert(rowI < Base::numRows);
    const auto finder = IJ2aI[rowI].find(colI);
    assert(finder != IJ2aI[rowI].end());
    Base::a(finder->second) = val;
    coefMtr.valuePtr()[finder->second] = val;
}

void EigenLibSolver::addCoeff(int rowI, int colI, double val){
    assert(rowI < Base::numRows);
    const auto finder = IJ2aI[rowI].find(colI);
    assert(finder != IJ2aI[rowI].end());
    Base::a(finder->second) += val;
    coefMtr.valuePtr()[finder->second] += val;
}

void EigenLibSolver::setZero(void){
    Base::setZero();
    memcpy(coefMtr.valuePtr(), Base::a.data(), Base::a.size()*sizeof(Base::a(0)));
}

void EigenLibSolver::setUnit_row(int rowI){
    assert(rowI < Base::numRows);
    for(const auto &item : Base::IJ2aI[rowI]){
        double tmp = (rowI == item.first);
        Base::a(item.second) = tmp;
        coefMtr.valuePtr()[item.second] = tmp;
    }
}

void EigenLibSolver::setZero_col(int colI, const std::set<int> &rowVIs){
    assert(colI < Base::numRows);
    for(const auto &rowI : rowVIs){
        assert(rowI < Base::numRows);
        const auto finder = Base::IJ2aI[rowI].find(colI);
        if(finder != Base::IJ2aI[rowI].end()){
            Base::a(finder->second) = 0.0;
            coefMtr.valuePtr()[finder->second] = 0.0;
        }
    }
}

} // namespace