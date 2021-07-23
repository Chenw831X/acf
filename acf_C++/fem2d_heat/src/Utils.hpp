//
// Utils.hpp
// acf_C++/fem2d_heat
//
// Created by Wei Chen on 7/23/21
//

#ifndef Utils_hpp
#define Utils_hpp

#include <Eigen/Eigen>

#include <string>

namespace fem2dHeat{

// a static class which implements some assistant functions
class Utils{
public:
    static bool readNodes(const std::string &filePath, Eigen::MatrixXd &V);
    static bool readTriangles(const std::string &filePath, Eigen::MatrixXi &F);
    static bool readDBC(const std::string &filePath, Eigen::MatrixXi &DBC);
    static bool readNBC(const std::string &filePath, Eigen::MatrixXi &NBC);
};

} // namespace

#endif // Utils_hpp