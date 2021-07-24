//
// LinSysSolver.hpp
// acf_C++/fem2d_heat
//
// Created by Wei Chen on 7/24/21
//

#ifndef LinSysSolver_hpp
#define LinSysSolver_hpp

#include <Eigen/Eigen>

#include <vector>
#include <map>
#include <set>

namespace fem2dHeat{

class LinSysSolver{
protected:
    int numRows, nnz;
    Eigen::VectorXi ia, ja;
    Eigen::VectorXd a;
    std::vector<std::map<int, int>> IJ2aI;

public:
    virtual ~LinSysSolver(void){}

public: // API
    virtual void set_pattern(const std::vector<std::set<int>> &vNeighbor){
        numRows = vNeighbor.size();
        nnz = 0;
        for(int vI=0; vI<numRows; ++vI){
            nnz += vNeighbor[vI].size();
        }

        ia.setZero(numRows+1);
        ia(0) = 0;
        ja.setZero(nnz);
        a.setZero(nnz);
        IJ2aI.resize(numRows);

        int cnt = 0;
        for(int vI=0; vI<numRows; ++vI){
            ia(vI+1) = ia(vI) + vNeighbor[vI].size();

            for(const auto &nbVI : vNeighbor[vI]){
                ja[cnt] = nbVI;
                IJ2aI[vI][nbVI] = cnt;
                ++cnt;
            }
        }
        assert(cnt == nnz);
    }

    virtual void analyze_pattern(void) = 0;

    virtual void factorize(void) = 0;

    virtual void solve(const Eigen::VectorXd &rhs, Eigen::VectorXd &res) = 0;

    virtual void setCoeff(int rowI, int colI, double val){
        assert(rowI < numRows);
        const auto finder = IJ2aI[rowI].find(colI);
        assert(finder != IJ2aI[rowI].end());
        a(finder->second) = val;
    }

    virtual void addCoeff(int rowI, int colI, double val){
        assert(rowI < numRows);
        const auto finder = IJ2aI[rowI].find(colI);
        assert(finder != IJ2aI[rowI].end());
        a(finder->second) += val;
    }

    virtual void setZero(void){
        a.setZero();
    }

    virtual void multiply(const Eigen::VectorXd &x, Eigen::VectorXd &Ax){
        assert(x.size() == numRows);
        Ax.setZero(numRows);
        for(int vI=0; vI<numRows; ++vI){
            for(const auto &item : IJ2aI[vI]){
                Ax(vI) += a[item.second] * x(item.first);
            }
        }
    }
    
    virtual void setUnit_row(int rowI){
        assert(rowI < numRows);
        for(const auto &item : IJ2aI[rowI]){
            a[item.second] = (rowI == item.first);
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

};

} // namespace

#endif // LinSysSolver_hpp