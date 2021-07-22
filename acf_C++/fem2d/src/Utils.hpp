//
// Utils.hpp
// acf_C++/fem2d
//
// Created by Wei Chen on 7/12/21
//

#ifndef Utils_hpp
#define Utils_hpp

#include <Eigen/Eigen>
#include <iostream>
#include <fstream>

namespace fem2d{

// a static class implementing some assitant functions
class Utils{
public:
    static bool readNodes(const std::string& filePath, Eigen::MatrixXd& V);
    static bool readTriangles(const std::string& filePath, Eigen::MatrixXi& F3);
    static bool readParallelograms(const std::string& filePath, Eigen::MatrixXi& F4);
    static bool readDirichlet(const std::string& filePath, Eigen::MatrixXi& DBC);
    static bool readNeumann(const std::string& filePath, Eigen::MatrixXi& NBC);
};

} // namespace

#endif //Utils_hpp