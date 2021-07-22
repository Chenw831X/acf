//
// LinSysSolver.hpp
// acf_C++/fem2d
//
// Created by Wei Chen on 7/10/21
//

#ifndef LinSysSolver_hpp
#define LinSysSolver_hpp

#include <Eigen/Eigen>
#include <iostream>
#include <set>
#include <map>
#include <vector>

namespace fem2d{

class LinSysSolver{
protected:
    int numRows, nnz;
    Eigen::VectorXi ia, ja;
    std::vector<std::map<int, int>> IJ2aI;
    Eigen::VectorXd a;

public:
    virtual ~LinSysSolver(void){};

    virtual void set_pattern(const std::vector<std::set<int>>& vNeighbor){
        numRows = static_cast<int>(vNeighbor.size());
        nnz = 0;
        for(int vI=0; vI<numRows; ++vI){
            nnz += static_cast<int>(vNeighbor[vI].size());
        }
        ia.resize(0);
        ja.resize(0);
        IJ2aI.resize(0);
        a.resize(0);

        ia.resize(numRows+1);
        ia(0) = 0;
        ja.resize(nnz);
        IJ2aI.resize(numRows);
        a.resize(nnz);

        for(int vI=0; vI<numRows; ++vI){
            int num = static_cast<int>(vNeighbor[vI].size());
            int cnt = 0;
            ia(vI+1) = ia(vI) + num;
            for(const auto& nbVI : vNeighbor[vI]){
                ja[ia(vI)+cnt] = nbVI;
                IJ2aI[vI][nbVI] = ia(vI) + cnt;
                ++cnt;
            }
        }

        // print_for_debug();
    }

    virtual void analyze_pattern(void) = 0;
    virtual void factorize(void) = 0;
    virtual void solve(const Eigen::VectorXd &rhs, Eigen::VectorXd &result) = 0;

    virtual void setCoeff(int rowI, int colI, double val){
        assert(rowI < IJ2aI.size());
        const auto finder = IJ2aI[rowI].find(colI);
        assert(finder != IJ2aI[rowI].end());
        a[finder->second] = val;
    }

    virtual void addCoeff(int rowI, int colI, double val){
        assert(rowI < IJ2aI.size());
        const auto finder = IJ2aI[rowI].find(colI);
        assert(finder != IJ2aI[rowI].end());
        a[finder->second] += val;
    }

    virtual void setZero(void){
        a.setZero();
    }

    virtual void multiply(const Eigen::VectorXd &x, Eigen::VectorXd &Ax){
        assert(x.size() == numRows);
        Ax.setZero(numRows);
        for(int rowI=0; rowI<numRows; ++rowI){
            for(const auto &colI : IJ2aI[rowI]){
                Ax(rowI) += a(colI.second) * x(colI.first);
            }
        }
    }

    virtual void setUnit_row(int rowI){
        assert(rowI < numRows);
        for(const auto &colIter : IJ2aI[rowI]){
            a[colIter.second] = (rowI==colIter.first);
        }
    }

    virtual void setZero_col(int colI, const std::set<int> &rowVIs){
        assert(colI < numRows);
        for(const auto &rowI : rowVIs){
            assert(rowI < numRows);
            const auto finder = IJ2aI[rowI].find(colI);
            if(finder != IJ2aI[rowI].end()){
                a[finder->second] = 0.0;
            }
        }
    }

    virtual void print_for_debug(void) const{
        std::cout << "ia:" << std::endl;
        std::cout << ia << std::endl;
        std::cout << "ja:" << std::endl;
        std::cout << ja << std::endl;
        std::cout << "a:" << std::endl;
        std::cout << a << std::endl;
        for(int vI=0; vI<numRows; ++vI){
            std::cout << vI << ":";
            for(const auto& item : IJ2aI[vI]){
                std::cout << " (" << item.first << "," << item.second << "),";
            }
            std::cout << std::endl;
        }
    }
};

} // namespace fem2d

#endif // LinSysSolver_hpp