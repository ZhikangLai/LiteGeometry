#include "PointTool.h"
//************************************************************************************************************************//
bool isPointInPolygon2D(
    const Eigen::MatrixX2d& polygonVertices,
    const Eigen::RowVector2d& point,
    bool checkBoundary,
    bool needClosePolygon
) {
    Eigen::MatrixX2d polygon;
    const Eigen::MatrixX2d* polygonPtr = nullptr;
    if (needClosePolygon) {
        polygon = generateClosedPolygon(polygonVertices);
        polygonPtr = &polygon;
    }
    else {
        polygonPtr = &polygonVertices;
    }

    Eigen::Index numVertices = polygonPtr->rows() - 1;

    bool inside = false;
    for (Eigen::Index i = 0; i < numVertices; ++i) {
        const auto& p1 = polygonPtr->row(i);
        const auto& p2 = polygonPtr->row(i + 1);

        if (isPointOnLine2D(Segment2D{ p1, p2 }, point)) {
            return checkBoundary;
        }

        const double y_min = std::min(p1.y(), p2.y());
        const double y_max = std::max(p1.y(), p2.y());
        if (point.y() < y_min || point.y() >= y_max) continue;

        const double lhs = (point.y() - p1.y()) * (p2.x() - p1.x());
        const double rhs = (point.x() - p1.x()) * (p2.y() - p1.y());
        if (lhs > rhs) inside = !inside;
    }
    return inside;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
bool isPointInPolygon3D(
    const Eigen::MatrixX3d& polygonVertices,
    const Eigen::RowVector3d& point,
    bool checkBoundary,
    bool needClosePolygon
) {
    if (!isCoplanar(polygonVertices)) return false;

    Eigen::MatrixX3d polygon;
    const Eigen::MatrixX3d* polygonPtr = nullptr;
    if (needClosePolygon) {
        polygon = generateClosedPolygon(polygonVertices);
        polygonPtr = &polygon;
    }
    else {
        polygonPtr = &polygonVertices;
    }

    const Eigen::RowVector3d& refPoint = polygonPtr->row(0);
    Eigen::RowVector3d faceNormal = computePlaneNormal(polygonVertices).first;
    const Eigen::RowVector3d deltaPoint = point - refPoint;
    double pointPlaneDist = std::abs(faceNormal.dot(deltaPoint));
    if (pointPlaneDist > epsilon) return false;

    Eigen::Matrix<double, 3, 2> T = compute3Dto2DTransformMatrix(faceNormal);

    Eigen::RowVector2d projectedPoint = deltaPoint * T;
    Eigen::MatrixX2d projectedVertices = (polygonPtr->rowwise() - refPoint) * T;

    return isPointInPolygon2D(projectedVertices, projectedPoint, checkBoundary, false);
}
//************************************************************************************************************************//

//************************************************************************************************************************//
bool isPointInPolyhedron(
    const Eigen::MatrixX3d& polyhedronVertices,
    const Eigen::RowVector3d& point,
    bool checkBoundary
) {
    Mesh polyhedron_mesh = buildConvexHullMesh(polyhedronVertices);
    auto face_iterators = faces(polyhedron_mesh);
    Tree tree(
        std::make_move_iterator(face_iterators.first),
        std::make_move_iterator(face_iterators.second),
        polyhedron_mesh);
    tree.build();
    tree.accelerate_distance_queries();
    Point_inside insideDetecter(std::move(tree));

    CGAL::Bounded_side result = insideDetecter(Point_3(point.x(), point.y(), point.z()));
    if (checkBoundary) {
        return result != CGAL::Bounded_side::ON_UNBOUNDED_SIDE;
    }
    else {
        return result == CGAL::Bounded_side::ON_BOUNDED_SIDE;
    }
}
//************************************************************************************************************************//