#ifndef DATAPREPARE_H
#define DATAPREPARE_H
#include "Utility.h"

struct EFData {
    Eigen::MatrixX3d polyhedronVertices = Eigen::MatrixX3d::Zero(8, 3);
    Eigen::Matrix<double, 8, 3> obbBoxVertices = Eigen::Matrix<double, 8, 3>::Zero();
    Eigen::RowVector3d obbBoxCenter = Eigen::RowVector3d::Zero();
    Eigen::Matrix<double, 8, 3> obbBoxExpand1Vertices = Eigen::Matrix<double, 8, 3>::Zero();
};

struct SCData {
    Eigen::MatrixX3d polyhedronVertices = Eigen::MatrixX3d::Zero(8, 3);
    Eigen::Matrix<double, 8, 3> obbBoxVertices = Eigen::Matrix<double, 8, 3>::Zero();
    Eigen::RowVector3d obbBoxCenter = Eigen::RowVector3d::Zero();
    Eigen::Matrix<double, 8, 3> obbBoxShrink1Vertices = Eigen::Matrix<double, 8, 3>::Zero();

};

std::pair<std::unordered_map<size_t, SCData>, std::unordered_map<size_t, EFData>> getDataSet(const std::string& filepath);
#endif 