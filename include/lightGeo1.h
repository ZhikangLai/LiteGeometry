#ifndef LIGHTGEO_H_
#define LIGHTGEO_H_
#include "Utility.h"
static constexpr std::array<int, 8> reorderIndexForCGAL = { 0,1,2,3,5,6,7,4 };

Mesh buildConvexHullMesh(const Eigen::MatrixX3d& polyhedronVertices);

Eigen::Matrix<double, 3, 2> computeProjectionMatrix(const Eigen::RowVector3d& planeNormal);

// Geometric tool : Point 
bool isCoplanar(const Eigen::MatrixX3d& points);

bool isPointOnLine2D(
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    const Eigen::RowVector2d& point,
    int lineType = 1
);

bool isPointOnLine3D(
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    const Eigen::RowVector3d& point,
    int lineType = 1
);

bool isPointInPolygon2D(
    const Eigen::MatrixX2d& polygonVertices,
    const Eigen::RowVector2d& point,
    bool checkBoundary = true,
    bool needClosePolygon = true
);

bool isPointInPolygon3D(
    const Eigen::MatrixX3d& polygonVertices,
    const Eigen::RowVector3d& point,
    bool checkBoundary = true,
    bool needClosePolygon = true
);

bool isPointInPolyhedron(
    const Eigen::MatrixX3d& _polyhedron,
    const Eigen::RowVector3d& point,
    bool checkBoundary = true
);
// Geometric tool : Point


// Geometric tool : Lines
bool isLinesIntersection2D(
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    const Eigen::RowVector2d& C,
    const Eigen::RowVector2d& D,
    Eigen::RowVector2d& intersection,
    int lineType = 1
);

bool isLinesIntersection2D(
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    const Eigen::RowVector2d& C,
    const Eigen::RowVector2d& D,
    int lineType = 1
);

std::tuple<Eigen::RowVector3d, Eigen::RowVector3d, double> computeLinesDistance(
    const Eigen::RowVector3d& A, const Eigen::RowVector3d& B,
    const Eigen::RowVector3d& C, const Eigen::RowVector3d& D,
    int lineType = 1
);

bool isLinesIntersection3D(
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    const Eigen::RowVector3d& C,
    const Eigen::RowVector3d& D,
    std::tuple<Eigen::RowVector3d, Eigen::RowVector3d, double>& output,
    double threshold = 0.0,
    int lineType = 1
);

bool isLinesIntersection3D(
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    const Eigen::RowVector3d& C,
    const Eigen::RowVector3d& D,
    double threshold = 0.0,
    int lineType = 1
);

bool isLinePolygonIntersection2D(
    const Eigen::MatrixX2d& polygonVertices,
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    std::vector<Eigen::RowVector2d>& intersections,
    int lineType = 1,
    bool needClosePolygon = true
);

bool isLinePolygonIntersection2D(
    const Eigen::MatrixX2d& polygonVertices,
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    int lineType = 1,
    bool needClosePolygon = true
);

bool isLinePolygonIntersection3D(
    const Eigen::MatrixX3d& polygonVertices,
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    std::vector<Eigen::RowVector3d>& intersections,
    int lineType,
    bool checkBoundary = true,
    bool needClosePolygon = true
);

bool isLinePolygonIntersection3D(
    const Eigen::MatrixX3d& polygonVertices,
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    int lineType = 1,
    bool checkBoundary = true,
    bool needClosePolygon = true
);


bool isLinePolyhedronIntersection(
    const Eigen::MatrixX3d& polyhedronVertices,
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    std::vector<Eigen::RowVector3d>& intersections,
    int lineType = 1
);

bool isLinePolyhedronIntersection(
    const Eigen::MatrixX3d& polyhedronVertices,
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    int lineType = 1
);

// Geometric tool : Lines

// Geometric tool : polygon and polyhedron
Eigen::MatrixXd generateClosedPolygon(const Eigen::MatrixXd& polygonVertices);

std::vector<Eigen::RowVector2d> filterPointsByPolygon(
    const Eigen::MatrixX2d& polygonVertices,
    const std::vector<Eigen::RowVector2d>& points,
    bool removeBoundary = true,
    bool removeInside = true,
    bool needClosePolygon = true
);

std::vector<Eigen::RowVector3d> filterPointsByPolyhedron(
    const Eigen::MatrixX3d& polyhedronVertices,
    const std::vector<Eigen::RowVector3d>& points,
    bool removeBoundary = true,
    bool removeInside = true
);


std::pair<Eigen::RowVector3d, Eigen::RowVector3d> computePlaneNormal(
    const Eigen::MatrixX3d& planePoints
);
// Geometric tool : polygon and polyhedron


// Minimum Bounding Shapes
std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computePCARect(
    const Eigen::MatrixX2d& points,
    bool useConvexHull = false
);

std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computePCABox(
    const Eigen::MatrixX3d& points,
    bool useConvexHull = false
);

std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computeMinBoundRect(
    const Eigen::MatrixX2d& points
);

std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computeMinBoundBox(
    const Eigen::MatrixX3d& points
);
// Minimum Bounding Shapes



// Projection & Coordinate Transformations
Eigen::MatrixX2d WorldToCameraImageCoords(
    const Eigen::MatrixX3d& targetVertices,
    const Eigen::RowVector3d& cameraDirPoint,
    const Eigen::RowVector3d& cameraPosition,
    double FL
);

Eigen::MatrixX2d WorldToCameraImageCoords(
    const Eigen::MatrixX3d& targetVertices,
    const Eigen::RowVector3d& cameraPosition,
    double FL
);
// Projection & Coordinate Transformations
#endif