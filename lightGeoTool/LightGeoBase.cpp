#include "LightGeoBase.h"

Mesh buildConvexHullMesh(
    const Eigen::MatrixX3d& polyhedronVertices
) {
    Eigen::Index numPoints = polyhedronVertices.rows();
    std::vector<Point_3> cgalPoints;
    cgalPoints.reserve(numPoints);
    for (Eigen::Index i = 0; i < numPoints; ++i) {
        cgalPoints.emplace_back(polyhedronVertices(i, 0), polyhedronVertices(i, 1), polyhedronVertices(i, 2));
    }
    Mesh mesh;
    CGAL::convex_hull_3(cgalPoints.begin(), cgalPoints.end(), mesh);
    return mesh;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
Eigen::MatrixXd generateClosedPolygon(
    const Eigen::MatrixXd& polygonVertices
) {
    Eigen::Index numRows = polygonVertices.rows();
    Eigen::Index numCols = polygonVertices.cols();
    if (numRows < 3) return polygonVertices;

    if (numCols != 2 && numCols != 3) {
        throw std::invalid_argument("Input must be 2D or 3D points");
    }
    if (numCols == 3 && !isCoplanar(polygonVertices)) {
        throw std::runtime_error("Non-coplanar 3D points detected");
    }

    Eigen::RowVectorXd center = polygonVertices.colwise().mean();
    std::vector<std::pair<Eigen::Index, double>> indexedAngles;
    indexedAngles.reserve(numRows);
    for (Eigen::Index i = 0; i < numRows; ++i) {
        Eigen::RowVectorXd vec = polygonVertices.row(i) - center;
        double angle = std::atan2(vec.y(), vec.x());
        if (angle < 0) angle += TWO_PI;
        indexedAngles.emplace_back(i, angle);
    }

    std::sort(indexedAngles.begin(), indexedAngles.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });

    Eigen::MatrixXd closedPolygon(numRows, numCols);
    for (Eigen::Index k = 0; k < numRows; ++k) {
        closedPolygon.row(k) = polygonVertices.row(indexedAngles[k].first);
    }

    if (closedPolygon.row(0).isApprox(closedPolygon.row(numRows - 1))) {
        return closedPolygon;
    }
    else {
        closedPolygon.conservativeResize(numRows + 1, Eigen::NoChange);
        closedPolygon.row(numRows) = closedPolygon.row(0);
        return closedPolygon;
    }

    return closedPolygon;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
std::pair<Eigen::RowVector3d, Eigen::RowVector3d> computePlaneNormal(
    const Eigen::MatrixX3d& planePoints
) {
    if (planePoints.rows() < 3) {
        throw std::invalid_argument("At least three points are required to define a plane.");
    }
    if (!isCoplanar(planePoints)) {
        throw std::invalid_argument("The provided points are not coplanar.");
    }
    Eigen::RowVector3d v1 = planePoints.row(1) - planePoints.row(0);
    Eigen::RowVector3d v2 = planePoints.row(2) - planePoints.row(1);
    Eigen::RowVector3d normal = v1.cross(v2).normalized();
    return { normal, -normal };
}
//************************************************************************************************************************//

//************************************************************************************************************************//
Eigen::Matrix<double, 3, 2> compute3Dto2DTransformMatrix(
    const Eigen::RowVector3d& planeNormal
) {
    Eigen::RowVector3d absNormal = planeNormal.cwiseAbs();
    int minIdx;
    absNormal.minCoeff(&minIdx);
    Eigen::RowVector3d v1;
    switch (minIdx) {
    case 0:
        v1 = Eigen::RowVector3d(0.0, -planeNormal.z(), planeNormal.y()).normalized();
        break;
    case 1:
        v1 = Eigen::RowVector3d(-planeNormal.z(), 0.0, planeNormal.x()).normalized();
        break;
    case 2:
        v1 = Eigen::RowVector3d(-planeNormal.y(), planeNormal.x(), 0.0).normalized();
        break;
    }

    Eigen::RowVector3d v2 = planeNormal.cross(v1);
    Eigen::Matrix<double, 3, 2> T;
    T << v1.transpose(), v2.transpose();
    return T;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
bool isCoplanar(const Eigen::MatrixX3d& points) {
    Eigen::MatrixX3d centered = points.bottomRows(points.rows() - 1).rowwise() - points.row(0);
    Eigen::JacobiSVD<Eigen::MatrixX3d> svd(centered, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Eigen::VectorXd& singularValues = svd.singularValues();
    return singularValues.size() < 3 || singularValues(2) <= epsilon;
}
//************************************************************************************************************************//