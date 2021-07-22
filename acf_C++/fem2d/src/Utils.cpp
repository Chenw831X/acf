//
// Utils.cpp
// acf_C++/fem2d
//
// Created by Wei Chen on 7/12/21
//

#include "Utils.hpp"

namespace fem2d{
    
bool Utils::readNodes(const std::string& filePath, Eigen::MatrixXd& V){
    FILE* in = fopen(filePath.c_str(), "r");
    if(!in){
        return false;
    }

    V.resize(0, 2);
    int bypass;
    double x, y;
    while(~fscanf(in, "%d%lf%lf", &bypass, &x, &y)){
        V.conservativeResize(V.rows()+1, 2);
        V(V.rows()-1, 0) = x;
        V(V.rows()-1, 1) = y;
    }

    fclose(in);
    return true;
}

bool Utils::readTriangles(const std::string& filePath, Eigen::MatrixXi& F3){
    FILE* in = fopen(filePath.c_str(), "r");
    if(!in){
        return false;
    }

    F3.resize(0, 3);
    int bypass, v1, v2, v3;
    while(~fscanf(in, "%d%d%d%d", &bypass, &v1, &v2, &v3)){
        F3.conservativeResize(F3.rows()+1, 3);
        F3(F3.rows()-1, 0) = v1-1;
        F3(F3.rows()-1, 1) = v2-1;
        F3(F3.rows()-1, 2) = v3-1;
    }

    fclose(in);
    return true;
}

bool Utils::readParallelograms(const std::string& filePath, Eigen::MatrixXi& F4){
    FILE* in = fopen(filePath.c_str(), "r");
    if(!in){
        return false;
    }

    F4.resize(0, 4);
    int bypass, v1, v2, v3, v4;
    while(~fscanf(in, "%d%d%d%d%d", &bypass, &v1, &v2, &v3, &v4)){
        F4.conservativeResize(F4.rows()+1, 4);
        F4(F4.rows()-1, 0) = v1-1;
        F4(F4.rows()-1, 1) = v2-1;
        F4(F4.rows()-1, 2) = v3-1;
        F4(F4.rows()-1, 3) = v4-1;
    }

    fclose(in);
    return true;
}

bool Utils::readDirichlet(const std::string& filePath, Eigen::MatrixXi& DBC){
    FILE* in = fopen(filePath.c_str(), "r");
    if(!in){
        return false;
    }

    DBC.resize(0, 2);
    int bypass, v1, v2;
    while(~fscanf(in, "%d%d%d", &bypass, &v1, &v2)){
        DBC.conservativeResize(DBC.rows()+1, 2);
        DBC(DBC.rows()-1, 0) = v1-1;
        DBC(DBC.rows()-1, 1) = v2-1;
    }

    fclose(in);
    return true;
}

bool Utils::readNeumann(const std::string& filePath, Eigen::MatrixXi& NBC){
    FILE* in = fopen(filePath.c_str(), "r");
    if(!in){
        return false;
    }

    NBC.resize(0, 2);
    int bypass, v1, v2;
    while(~fscanf(in, "%d%d%d", &bypass, &v1, &v2)){
        NBC.conservativeResize(NBC.rows()+1, 2);
        NBC(NBC.rows()-1, 0) = v1-1;
        NBC(NBC.rows()-1, 1) = v2-1;
    }

    fclose(in);
    return true;
}

} // namespace