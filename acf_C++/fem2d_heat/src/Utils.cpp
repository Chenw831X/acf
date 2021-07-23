//
// Utils.cpp
// acf_C++/fem2d_heat
//
// Created by Wei Chen on 7/23/21
//

#include "Utils.hpp"

namespace fem2dHeat{

bool Utils::readNodes(const std::string &filePath, Eigen::MatrixXd &V){
    FILE *in = fopen(filePath.c_str(), "r");
    if(!in){
        return false;
    }

    V.resize(0, 2);
    double bypass, x, y;
    while(~fscanf(in, "%le%le%le", &bypass, &x, &y)){
        V.conservativeResize(V.rows()+1, 2);
        V(V.rows()-1, 0) = x;
        V(V.rows()-1, 1) = y;
    }

    fclose(in);
    return true;
}

bool Utils::readTriangles(const std::string &filePath, Eigen::MatrixXi &F){
    FILE* in = fopen(filePath.c_str(), "r");
    if(!in){
        return false;
    }

    F.resize(0, 3);
    double bypass, v1, v2, v3;
    while(~fscanf(in, "%le%le%le%le", &bypass, &v1, &v2, &v3)){
        F.conservativeResize(F.rows()+1, 3);
        F(F.rows()-1, 0) = (int)v1;
        F(F.rows()-1, 1) = (int)v2;
        F(F.rows()-1, 2) = (int)v3;
    }

    fclose(in);
    return true;
}

bool Utils::readDBC(const std::string &filePath, Eigen::MatrixXi &DBC){
    FILE *in = fopen(filePath.c_str(), "r");
    if(!in){
        return false;
    }

    DBC.resize(0, 2);
    double bypass, v1, v2;
    while(~fscanf(in, "%le%le%le", &bypass, &v1, &v2)){
        DBC.conservativeResize(DBC.rows()+1, 2);
        DBC(DBC.rows()-1, 0) = (int)v1;
        DBC(DBC.rows()-1, 1) = (int)v2;
    }

    fclose(in);
    return true;
}

bool Utils::readNBC(const std::string &filePath, Eigen::MatrixXi &NBC){
    FILE *in = fopen(filePath.c_str(), "r");
    if(!in){
        return false;
    }

    NBC.resize(0, 2);
    double bypass, v1, v2;
    while(~fscanf(in, "%le%le%le", &bypass, &v1, &v2)){
        NBC.conservativeResize(NBC.rows()+1, 2);
        NBC(NBC.rows()-1, 0) = (int)v1;
        NBC(NBC.rows()-1, 1) = (int)v2;
    } 

    fclose(in);
    return true;
}

} // namespace