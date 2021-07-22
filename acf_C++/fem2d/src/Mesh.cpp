//
// Mesh.cpp
// acf_C++/fem2d
//
// Created by Wei Chen on 7/12/21
//

#include "Mesh.hpp"
#include "Utils.hpp"

#include <Eigen/Eigen>
#include <iostream>

namespace fem2d{

Mesh::Mesh(void){
    if(!Utils::readNodes("../input/coordinates.dat", V)){
        std::cout << "Logging: [fem2d] error on reading coordinates" << std::endl;
        exit(-1);
    }
    if(!Utils::readTriangles("../input/elements3.dat", F3)){
        std::cout << "Logging: [fem2d] error on reading triangles" << std::endl;
        exit(-1);
    }
    if(!Utils::readParallelograms("../input/elements4.dat", F4)){
        std::cout << "Logging: [fem2d] error on reading parallelograms" << std::endl;
        exit(-1);
    }
    if(!Utils::readDirichlet("../input/dirichlet.dat", DBC)){
        std::cout << "Logging: [fem2d] error on reading dirichlet boundary conditions" << std::endl;
        exit(-1);
    }
    if(!Utils::readNeumann("../input/neumann.dat", NBC)){
        std::cout << "Logging: [fem2d] error on reading neumann boundary conditions" << std::endl;
        exit(-1);
    }
    std::cout << "Logging: [fem2d] read mesh completed!" << std::endl;

    computeFeatures();
    std::cout << "Logging: [fem2d] compute features" << std::endl;

    linSysSolver = new EigenLibSolver();
    linSysSolver->set_pattern(vNeighbor);
    std::cout << "Logging: [fem2d] set pattern for linear solver" << std::endl;
    linSysSolver->analyze_pattern();
    std::cout << "Logging: [fem2d] analyze pattern" << std::endl;
    // print_for_debug();
}

Mesh::~Mesh(void){
    delete linSysSolver;
}

void Mesh::print_for_debug(void) const{
    /*
    std::cout << "coordinates:" << std::endl;
    std::cout << V << std::endl << std::endl;
    std::cout << "triangles:" << std::endl;
    std::cout << F3 << std::endl << std::endl;
    std::cout << "parallelograms:" << std::endl;
    std::cout << F4 << std::endl << std::endl;
    std::cout << "dirichlet:" << std::endl;
    std::cout << DBC << std::endl << std::endl;
    std::cout << "neumann:" << std::endl;
    std::cout << NBC << std::endl << std::endl;
    std::cout << "vNeighbor:" << std::endl;
    for(int vI=0; vI<V.rows(); ++vI){
        std::cout << vI << ":";
        for(const auto& nbVI : vNeighbor[vI]){
            std::cout << " " << nbVI;
        }
        std::cout << std::endl;
    }
    std::cout << "DBCV:" << std::endl;
    for(const auto& DBCI : DBCV){
        std::cout << DBCI << std::endl;
    }
    std::cout << "FreeV:" << std::endl;
    for(const auto &FreeI : FreeV){
        std::cout << FreeI << std::endl;
    }
    std::cout << "b" << std::endl;
    std::cout << b << std::endl;
    std::cout << "u=" << std::endl;
    std::cout << u << std::endl;
    */
}

bool Mesh::writeResult(void) const{
    FILE *out = fopen("../output/u.txt", "w");
    if(!out){
        return false;
    }

    for(int vI=0; vI<u.size(); ++vI)
        fprintf(out, "%f\n", u(vI));    

    fclose(out);
    return true;
}

void Mesh::computeFeatures(void){ // compute vNeighbor and DBCV
    // vNeighbor
    vNeighbor.resize(0);
    vNeighbor.resize(V.rows());

    for(int vI=0; vI<V.rows(); ++vI){
        vNeighbor[vI].insert(vI);
    }

    for(int triI=0; triI<F3.rows(); ++triI){
        const Eigen::Matrix<int, 1, 3>& triVInd = F3.row(triI);
        for(int vI=0; vI<3; ++vI){
            for(int vJ=vI+1; vJ<3; ++vJ){
                vNeighbor[triVInd(vI)].insert(triVInd(vJ));
                vNeighbor[triVInd(vJ)].insert(triVInd(vI));
            }
        }
    }

    for(int paraI=0; paraI<F4.rows(); ++paraI){
        const Eigen::Matrix<int, 1, 4>& paraVInd = F4.row(paraI);
        for(int vI=0; vI<4; ++vI){
            for(int vJ=vI+1; vJ<4; ++vJ){
                vNeighbor[paraVInd(vI)].insert(paraVInd(vJ));
                vNeighbor[paraVInd(vJ)].insert(paraVInd(vI));
            }
        }
    }

    // DBCV
    for(int DBCI=0; DBCI<DBC.rows(); ++DBCI){
        const Eigen::Matrix<int, 1, 2> &DBCVInd = DBC.row(DBCI);
        DBCV.insert(DBCVInd(0));
        DBCV.insert(DBCVInd(1));
    }

    // FreeV
    for(int vI=0; vI<V.rows(); ++vI){
        if(DBCV.find(vI) == DBCV.end()){
            FreeV.insert(vI);
        }
    }
}

void Mesh::assembly_stiffnessMat(void){ // assebly stiffness matrix
    linSysSolver->setZero();
    // assembly triangle elements' stiffness matrix
    for(int triI=0; triI<F3.rows(); ++triI){
        const Eigen::Matrix<int, 1, 3>& triVInd = F3.row(triI);
        Eigen::Matrix3d local;
        stima3(triI, local);
        for(int i=0; i<3; ++i){
            int rowI = triVInd(i);
            for(int j=0; j<3; ++j){
                int colI = triVInd(j);
                linSysSolver->addCoeff(rowI, colI, local(i, j));
            }
        }
    }
    // assembly parallelograms elements' stiffness matrix
    for(int paraI=0; paraI<F4.rows(); ++paraI){
        const Eigen::Matrix<int, 1, 4>& paraVInd = F4.row(paraI);
        Eigen::Matrix4d local;
        stima4(paraI, local);
        for(int i=0; i<4; ++i){
            int rowI = paraVInd(i);
            for(int j=0; j<4; ++j){
                int colI = paraVInd(j);
                linSysSolver->addCoeff(rowI, colI, local(i, j));
            }
        }
    }
}

void Mesh::assembly_rhs(void){ // assembly right hand side
    b.resize(V.rows());
    b.setZero();
    
    // Volume Forces
    for(int triI=0; triI<F3.rows(); ++triI){ 
        const Eigen::Matrix<int, 1, 3>& triVInd = F3.row(triI);

        Eigen::Matrix3d temp1;
        temp1.setZero();
        temp1(0, 0) = 1.0;
        temp1(1, 0) = V(triVInd(0), 0);
        temp1(2, 0) = V(triVInd(0), 1);
        temp1(0, 1) = 1.0;
        temp1(1, 1) = V(triVInd(1), 0);
        temp1(2, 1) = V(triVInd(1), 1);
        temp1(0, 2) = 1.0;
        temp1(1, 2) = V(triVInd(2), 0);
        temp1(2, 2) = V(triVInd(2), 1);

        Eigen::Matrix<double, 1, 2> Vmed;
        Vmed.setZero();
        Vmed += V.row(triVInd(0));
        Vmed += V.row(triVInd(1));
        Vmed += V.row(triVInd(2));
        Vmed /= 3.0;

        double res = temp1.determinant() * f(Vmed) / 6.0;
        b(triVInd(0)) += res;
        b(triVInd(1)) += res;
        b(triVInd(2)) += res;
    }

    for(int paraI=0; paraI<F4.rows(); ++paraI){
        const Eigen::Matrix<int, 1, 4> &paraVInd = F4.row(paraI);

        Eigen::Matrix3d temp1;
        temp1.setZero();
        temp1(0, 0) = 1.0;
        temp1(1, 0) = V(paraVInd(0), 0);
        temp1(2, 0) = V(paraVInd(0), 1);
        temp1(0, 1) = 1.0;
        temp1(1, 1) = V(paraVInd(1), 0);
        temp1(2, 1) = V(paraVInd(1), 1);
        temp1(0, 2) = 1.0;
        temp1(1, 2) = V(paraVInd(2), 0);
        temp1(2, 2) = V(paraVInd(2), 1);

        Eigen::Matrix<double, 1, 2> Vmed;
        Vmed.setZero();
        Vmed += V.row(paraVInd(0));
        Vmed += V.row(paraVInd(1));
        Vmed += V.row(paraVInd(2));
        Vmed += V.row(paraVInd(3));
        Vmed /= 4.0;

        double res = temp1.determinant() * f(Vmed) / 4.0;
        b(paraVInd(0)) += res;
        b(paraVInd(1)) += res;
        b(paraVInd(2)) += res;
        b(paraVInd(3)) += res;
    }

    // Neumann boundary conditions
    for(int NBCI=0; NBCI<NBC.rows(); ++NBCI){
        const Eigen::Matrix<int, 1, 2> &NBCVInd = NBC.row(NBCI);
        double length = (V.row(NBCVInd(0)) - V.row(NBCVInd(1))).norm();

        Eigen::Matrix<double, 1, 2> Vmed;
        Vmed.setZero();
        Vmed += V.row(NBCVInd(0));
        Vmed += V.row(NBCVInd(1));
        Vmed /= 2.0;

        double res = length * g(Vmed) / 2.0;
    }

    // Dirichlet boundary conditions
    u.resize(V.rows());
    u.setZero();
    u_d();    
    Eigen::VectorXd Ax;
    linSysSolver->multiply(u, Ax);
    b = b - Ax;
}

void Mesh::adjustByDirichlet(void){ // adjust A and b based on dirichlet boundary conditions
    for(const auto& DBCI : DBCV){ // adjust b
        b(DBCI) = u(DBCI);
    }

    for(const auto& DBCI : DBCV){ // adjust A
        linSysSolver->setZero_col(DBCI, FreeV);
    }
    for(const auto& DBCI : DBCV){
        linSysSolver->setUnit_row(DBCI);
    }

    // print_for_debug();
    // linSysSolver->print_for_debug();
}

void Mesh::solve(void){ // computation of the solution
    linSysSolver->factorize();
    linSysSolver->solve(b, u); 
    std::cout << "Logging: [fem2d] solve linear system" << std::endl;

    if(!writeResult()){
        std::cout << "Logging: [fem2d] error on writing u" << std::endl;
        exit(-1);
    }
    std::cout << "Logging: [fem2d] writing u to output/u.txt" << std::endl;
}

// <<< interfaces
void Mesh::stima3(int triI, Eigen::Matrix3d &local){
    const Eigen::Matrix<int, 1, 3>& triVInd = F3.row(triI);

    Eigen::Matrix3d temp1;
    temp1.setZero();
    temp1(0, 0) = 1.0;
    temp1(1, 0) = V(triVInd(0), 0);
    temp1(2, 0) = V(triVInd(0), 1);
    temp1(0, 1) = 1.0;
    temp1(1, 1) = V(triVInd(1), 0);
    temp1(2, 1) = V(triVInd(1), 1);
    temp1(0, 2) = 1.0;
    temp1(1, 2) = V(triVInd(2), 0);
    temp1(2, 2) = V(triVInd(2), 1);

    Eigen::Matrix<double, 3, 2> temp2;
    temp2.setZero();
    temp2(1, 0) = 1.0;
    temp2(2, 1) = 1.0;

    Eigen::Matrix<double, 3, 2> G = temp1.inverse() * temp2;
    local = temp1.determinant() / 2.0 * G * G.transpose();
}

void Mesh::stima4(int paraI, Eigen::Matrix4d &local){
    const Eigen::Matrix<int, 1, 4>& paraVInd = F4.row(paraI);

    Eigen::Matrix2d D_Phi;
    D_Phi.setZero();
    D_Phi(0, 0) = V(paraVInd(1), 0) - V(paraVInd(0), 0);
    D_Phi(1, 0) = V(paraVInd(1), 1) - V(paraVInd(0), 1);
    D_Phi(0, 1) = V(paraVInd(3), 0) - V(paraVInd(0), 0);
    D_Phi(1, 1) = V(paraVInd(3), 1) - V(paraVInd(0), 1);
    Eigen::Matrix2d B = (D_Phi.transpose() * D_Phi).inverse();
    double a = B(0, 0);
    double b = B(0, 1);
    double c = B(1, 1);
    double C1_00 = 3.0 * b + 2.0 * (a + c);
    double C1_01 = -2.0 * a + c;
    double C1_10 = C1_01;
    double C1_11 = -3.0 * b + 2.0 * (a + c);
    double C2_00 = -3.0 * b - (a + c);
    double C2_01 = a - 2.0 * c;
    double C2_10 = C2_01;
    double C2_11 = 3.0 * b - (a + c);
    local(0, 0) = C1_00;
    local(0, 1) = C1_01;
    local(1, 0) = C1_10;
    local(1, 1) = C1_11;
    local(0, 2) = C2_00;
    local(0, 3) = C2_01;
    local(1, 2) = C2_10;
    local(1, 3) = C2_11;
    local(2, 0) = C2_00;
    local(2, 1) = C2_01;
    local(3, 0) = C2_10;
    local(3, 1) = C2_11;
    local(2, 2) = C1_00;
    local(2, 3) = C1_01;
    local(3, 2) = C1_10;
    local(3, 3) = C1_11;
    local *= D_Phi.determinant() / 6.0;
}

double Mesh::f(const Eigen::Matrix<double, 1, 2> &Vmed){
    return 1.0;
}

double Mesh::g(const Eigen::Matrix<double, 1, 2> &Vmed){
    return 0.0;
}

void Mesh::u_d(void){
    for(const auto& DBCI : DBCV){
        u(DBCI) = 0.0;
    }
}

// >>> interfaces

} // namespace