#ifndef LINETOOL_H_
#define LINETOOL_H_
#include "LineTraits.h"
#include "LightGeoBase.h"

//************************************************************************************************************************//
/**
 * @brief Computes the closest points and shortest distance between two 3D linear primitives (segments, rays, or lines).
 * 
 * The function works with `Segment3D`, `Ray3D`, or `Line3D` types. These types must have members `P1` and `P2` of 
 * type `Eigen::RowVector3d` representing the start and end points of the primitive.
 * 
 * @tparam LineT: The type of the 3D linear primitives (e.g., `Segment3D`, `Ray3D`, `Line3D`).
 *               Must define members P1 and P2 (both Eigen::RowVector3d).
 * @param[in] geoAB: The first linear primitives.
 * @param[in] geoCD: The second linear primitives.
 * @return A tuple consisting of:
 *         - The closest point on `geoAB` (`Eigen::RowVector3d`)
 *         - The closest point on `geoCD` (`Eigen::RowVector3d`)
 *         - The shortest distance between the primitives (`double`)
 */
template<typename LineT>
std::tuple<Eigen::RowVector3d, Eigen::RowVector3d, double> computeLinesDistance(
    const LineT& geoAB, const LineT& geoCD)
{
    const Eigen::RowVector3d dirAB = geoAB.P2 - geoAB.P1;
    const Eigen::RowVector3d dirCD = geoCD.P2 - geoCD.P1;
    const Eigen::RowVector3d dirAC = geoCD.P1 - geoAB.P1;

    const double lenSqAB = dirAB.squaredNorm();
    const double lenSqCD = dirCD.squaredNorm();
    double t = 0.0, u = 0.0;
    if (lenSqAB < highEpsilon || lenSqCD < highEpsilon) {
        if (lenSqAB > highEpsilon) {
            t = dirAB.dot(dirAC) / lenSqAB;
        }
        else if (lenSqCD > highEpsilon) {
            u = -dirCD.dot(dirAC) / lenSqCD;
        }
    }
    else {
        const double AB_dot_CD = dirAB.dot(dirCD);
        const double denom = lenSqAB * lenSqCD - AB_dot_CD * AB_dot_CD;
        if (std::abs(denom) < highEpsilon) {// 平行或共线
            u = -dirCD.dot(dirAC) / lenSqCD;
            t = (u * AB_dot_CD + dirAB.dot(dirAC)) / lenSqAB;
        }
        else {// 一般情况
            t = (dirAB.dot(dirAC) * lenSqCD - dirCD.dot(dirAC) * AB_dot_CD) / denom;
            u = (t * AB_dot_CD - dirCD.dot(dirAC)) / lenSqCD;
        }
    }

    line_traits<LineT>::clampParams(t, u);

    const Eigen::RowVector3d closestOnAB = geoAB.P1 + dirAB * t;
    const Eigen::RowVector3d closestOnCD = geoCD.P1 + dirCD * u;

    return { closestOnAB, closestOnCD, (closestOnAB - closestOnCD).norm() };
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if two 2D linear primitives intersect and calculates the intersection point.
 * 
 * Supports `Segment2D`, `Ray2D`, or `Line2D` types with `P1` and `P2` members of type `Eigen::RowVector2d`.
 * 
 * @tparam LineT: The type of the 2D linear primitive; it must provide `Eigen::RowVector2d` members `P1` and `P2`.
 * @param[in] geoAB: The first linear primitive.
 * @param[in] geoCD: The second linear primitive.
 * @param[out] intersection: If the lines intersect, this parameter will hold the intersection point of type Eigen::RowVector2d.
 * @return Return `true` if an intersection is found, `false` otherwise.
 */
template<typename LineT>
bool isLinesIntersection2D(
    const LineT& geoAB,
    const LineT& geoCD,
    Eigen::RowVector2d& intersection
) {
    const Eigen::RowVector2d dirAB = geoAB.P2 - geoAB.P1;
    const Eigen::RowVector2d dirCD = geoCD.P2 - geoCD.P1;
    const Eigen::RowVector2d dirAC = geoCD.P1 - geoAB.P1;

    const double denom = dirAB.cross(dirCD);

    if (std::abs(denom) < epsilon) return false;

    const double tScaled = dirAC.cross(dirCD);
    const double uScaled = dirAC.cross(dirAB);

    if (!line_traits<LineT>::validIntersect(tScaled, uScaled, denom)) return false;

    intersection = geoAB.P1 + dirAB * (tScaled / denom);
    return true;
}

/**
 * @brief Checks whether two 2D linear primitives (segments, rays, or lines) intersect.
 *
 * The function works with `Segment2D`, `Ray2D`, or `Line2D` types, all of which must have `P1` and `P2` as members of type
 * Eigen::RowVector2d representing the start and end points of the linear primitives.
 *
 * @tparam LineT: The type of the 2D linear primitive; it must provide `Eigen::RowVector2d` members `P1` and `P2`.
 * @param[in] geoAB: The first linear primitive.
 * @param[in] geoCD: The second linear primitive.
 * @return Returns `true` if intersection exists; `false` otherwise.
 */
template<typename LineT>
bool isLinesIntersection2D(
    const LineT& geoAB,
    const LineT& geoCD
) {
    const Eigen::RowVector2d dirAB = geoAB.P2 - geoAB.P1;
    const Eigen::RowVector2d dirCD = geoCD.P2 - geoCD.P1;
    const Eigen::RowVector2d dirAC = geoCD.P1 - geoAB.P1;

    const double denom = dirAB.cross(dirCD);
    if (std::abs(denom) < epsilon) return false;

    const double tScaled = dirAC.cross(dirCD);
    const double uScaled = dirAC.cross(dirAB);

    if (line_traits<LineT>::validIntersect(tScaled, uScaled, denom)) {
        return true;
    }
    else {
        return false;
    }
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if two 3D linear primitives (segments, rays, or lines) intersect based on a distance threshold.
 *
 * The function works with `Segment3D`, `Ray3D`, or `Line3D` types, all of which must have `P1` and `P2` as members of type
 * Eigen::RowVector3d representing the start and end points of the linear primitives. It computes the shortest distance between 
 * two 3D linear primitives (e.g., Segment3D, Ray3D, Line3D) and checks if the distance is within a specified threshold.
 *
 * @tparam T1: The type of the first 3D linear primitive; it must provide `Eigen::RowVector3d` members `P1` and `P2`.
 * @tparam T2: The type of the second 3D linear primitive; it must provide `Eigen::RowVector3d` members `P1` and `P2`.
 * @param[in] geoAB: The first linear primitive.
 * @param[in] geoCD: The second linear primitive.
 * @param[out] output: If the shortest distance is less than or equal to the threshold, this parameter will hold
 * the closest points and the distance as a tuple (Eigen::RowVector3d, Eigen::RowVector3d, double).
 * @param[in] threshold: The maximum allowable distance for the primitives to be considered intersecting.
 * @return Returns `true` if distance ≤ threshold; `false` otherwise.
 */
template<typename T1, typename T2>
bool isLinesIntersection3D(
    const T1& geoAB,
    const T2& geoCD,
    std::tuple<Eigen::RowVector3d, Eigen::RowVector3d, double>& output,
    double threshold
) {
    output = computeLinesDistance(geoAB, geoCD);
    return std::get<2>(output) <= threshold + epsilon;
}

/**
 * @brief Checks if two 3D linear primitives (segments, rays, or lines) intersect based on a distance threshold.
 *
 * The function works with `Segment3D`, `Ray3D`, or `Line3D` types, all of which must have `P1` and `P2` as members of type
 * Eigen::RowVector3d representing the start and end points of the linear primitives. It computes the shortest distance between
 * two 3D linear primitives (e.g., Segment3D, Ray3D, Line3D) and checks if the distance is within a specified threshold.
 *
 * @tparam T1: The type of the first 3D linear primitive; it must provide `Eigen::RowVector3d` members `P1` and `P2`.
 * @tparam T2: The type of the second 3D linear primitive; it must provide `Eigen::RowVector3d` members `P1` and `P2`.
 * @param[in] geoAB: The first linear primitive.
 * @param[in] geoCD: The second linear primitive.
 * @param[in] threshold: The maximum allowable distance for the primitives to be considered intersecting.
 * @return Returns `true` if distance ≤ threshold; `false` otherwise.
 */
template<typename T1, typename T2>
bool isLinesIntersection3D(
    const T1& geoAB,
    const T2& geoCD,
    double threshold
) {
    return std::get<2>(computeLinesDistance(geoAB, geoCD)) <= threshold + epsilon;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if a 2D linear primitives (segments, rays, or lines) intersects with a polygon and calculates intersection points.
 *
 * The function works with `Segment2D`, `Ray2D`, or `Line2D` types, all of which must have `P1` and `P2` as members of type
 * Eigen::RowVector2d representing the start and end points of the linear primitives.
 *
 * @tparam LineT: The type of the 2D linear primitive; it must provide `Eigen::RowVector2d` members `P1` and `P2`.
 * @param[in] polygonVertices: A matrix of size `N×2` where each row represents a 2D vertex of the 2D polygon.
 *                             The vertices may not be ordered or closed. The function will automatically close and sort the polygon if necessary.
 * @param[in] geoAB: The 2D linear primitive  to check for intersection with the polygon.
 * @param[out] intersections: If the lines intersect the 2D polygon, this parameter will hold the intersection points, each of type `Eigen::RowVector2d`.
 * @param[in] needClosePolygon: If `true` (default), the function will automatically close the polygon by sorting its vertices
 *                             counter-clockwise and ensuring the polygon is closed.
 *                             If `false`, you must ensure that the `polygonVertices` is sorted and closed.
 *                             Setting this to `false` avoids the overhead of sorting and closing, improving performance when the polygon is already sorted and closed.
 * @return Returns `true` if at least one intersection is found; `false` otherwise.
 */
template<typename LineT>
bool isLinePolygonIntersection2D(
    const Eigen::MatrixX2d& polygonVertices,
    const LineT& geoAB,
    std::vector<Eigen::RowVector2d>& intersections,
    bool needClosePolygon = true
) {
    Eigen::MatrixX2d polygon;
    const Eigen::MatrixX2d* polygonPtr;
    if (needClosePolygon) {
        polygon = generateClosedPolygon(polygonVertices);
        polygonPtr = &polygon;
    }
    else {
        polygonPtr = &polygonVertices;
    }

    const Eigen::Index numVertices = polygonPtr->rows() - 1;
    intersections.reserve(2);
    for (Eigen::Index i = 0; i < numVertices; ++i) {
        const auto& p1 = polygonPtr->row(i);
        const auto& p2 = polygonPtr->row(i + 1);
        Eigen::RowVector2d intersection;
        if (isLinesIntersection2D(geoAB, LineT{ p1, p2 }, intersection)) {
            intersections.emplace_back(intersection);
        }
    }
    return !intersections.empty();
}

/**
 * @brief Checks if a 2D linear primitives (segments, rays, or lines) intersects with a 2D polygon.
 *
 * The function works with `Segment2D`, `Ray2D`, or `Line2D` types, all of which must have `P1` and `P2` as members of type
 * Eigen::RowVector2d representing the start and end points of the linear primitives. 
 *
 * @tparam LineT: The type of the 2D linear primitive; it must provide `Eigen::RowVector2d` members `P1` and `P2`.
 * @param[in] polygonVertices: A matrix of size `N×2` where each row represents a 2D vertex of the polygon.
 *                             The vertices may not be ordered or closed. The function will automatically close and sort the polygon if necessary.
 * @param[in] geoAB: The 2D linear primitive to check for intersection with the polygon.
 * @param[in] needClosePolygon: If `true` (default), the function will automatically close the 2D polygon by sorting its vertices
 *                             counter-clockwise and ensuring the polygon is closed.
 *                             If `false`, you must ensure that the `polygonVertices` is sorted and closed.
 *                             Setting this to `false` avoids the overhead of sorting and closing, improving performance when the polygon is already sorted and closed.
 * @return Returns `true` if at least one intersection is found; `false` otherwise.
 */
template<typename LineT>
bool isLinePolygonIntersection2D(
    const Eigen::MatrixX2d& polygonVertices,
    const LineT& geoAB,
    bool needClosePolygon = true
) {
    Eigen::MatrixX2d polygon;
    const Eigen::MatrixX2d* polygonPtr;
    if (needClosePolygon) {
        polygon = generateClosedPolygon(polygonVertices);
        polygonPtr = &polygon;
    }
    else {
        polygonPtr = &polygonVertices;
    }

    const Eigen::Index numVertices = polygonPtr->rows() - 1;
    for (Eigen::Index i = 0; i < numVertices; ++i) {
        const auto& p1 = polygonPtr->row(i);
        const auto& p2 = polygonPtr->row(i + 1);
        if (isLinesIntersection2D(geoAB, LineT{ p1, p2 })) {
            return true;
        }
    }
    return false;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if a 3D linear primitives (segments, rays, or lines) intersects with a 3D polygon and calculates intersection points.
 *
 * The function works with `Segment3D`, `Ray3D`, or `Line3D` types, all of which must have `P1` and `P2` as members of type
 * Eigen::RowVector3d representing the start and end points of the linear primitives.
 *
 * @tparam LineT: The type of the 3D linear primitive; it must provide `Eigen::RowVector3d` members `P1` and `P2`.
 * @param[in] polygonVertices: A matrix of size `N×3` where each row represents a 3D vertex of the 3D polygon.
 *                             The vertices may not be ordered or closed. The function will automatically close and sort the polygon if necessary.
 * @param[in] geoAB: The 3D linear primitive to check for intersection with the polygon.
 * @param[out] intersections: If the lines intersect the 3D polygon, this parameter will hold the intersection points, each of type `Eigen::RowVector3d`.
 * @param[in] checkBoundary: If `true` (default), the function considers intersections on the polygon's boundaries. 
 *                            If `false`, only intersections strictly inside the polygon are considered.
 * @param[in] needClosePolygon: If `true` (default), the function will automatically close the polygon by sorting its vertices
 *                             counter-clockwise and ensuring the polygon is closed.
 *                             If `false`, you must ensure that the `polygonVertices` is sorted and closed.
 *                             Setting this to `false` avoids the overhead of sorting and closing, improving performance when the polygon is already sorted and closed.
 * @return Returns `true` if at least one intersection is found; `false` otherwise.
 */
template<typename LineT>
bool isLinePolygonIntersection3D(
    const Eigen::MatrixX3d& polygonVertices,
    const LineT& geoAB,
    std::vector<Eigen::RowVector3d>& intersections,
    bool checkBoundary = true,
    bool needClosePolygon = true
) {
    if (!isCoplanar(polygonVertices)) return false;

    Eigen::MatrixX3d polygon;
    const Eigen::MatrixX3d* polygonPtr;
    if (needClosePolygon) {
        polygon = generateClosedPolygon(polygonVertices);
        polygonPtr = &polygon;
    }
    else {
        polygonPtr = &polygonVertices;
    }

    const Eigen::RowVector3d& refPoint = polygonPtr->row(0);
    Eigen::RowVector3d faceNormal = computePlaneNormal(*polygonPtr).first;
    const Eigen::RowVector3d dirAB = geoAB.P2 - geoAB.P1;
    const double denom = faceNormal.dot(dirAB);
    const Eigen::RowVector3d deltaA = geoAB.P1 - refPoint;
    const bool isParallel = std::abs(denom) < epsilon;
    if (checkBoundary && isParallel) {
        const Eigen::RowVector3d deltaB = geoAB.P2 - refPoint;
        if (std::abs(faceNormal.dot(deltaA)) > epsilon ||
            std::abs(faceNormal.dot(deltaB)) > epsilon) return false;

        Eigen::Matrix<double, 3, 2> T = compute3Dto2DTransformMatrix(faceNormal);
        Eigen::MatrixX3d deltaVertices = polygonPtr->rowwise() - refPoint;
        std::vector<Eigen::RowVector2d> intersections2D;

        using Line2DType = typename line_traits<LineT>::Line2DType;
        if (isLinePolygonIntersection2D(deltaVertices * T, Line2DType{ deltaA * T, deltaB * T }, intersections2D, false)) {
            for (const auto& p2d : intersections2D) {
                intersections.emplace_back(refPoint + p2d * T.transpose());
            }
            return true;
        }
        else {
            return false;
        }
    }

    if (isParallel) return false;

    const double t = faceNormal.dot(-deltaA) / denom;
    if (!line_traits<LineT>::validT(t)) return false;

    Eigen::RowVector3d intersection = geoAB.P1 + t * dirAB;
    Eigen::Matrix<double, 3, 2> T = compute3Dto2DTransformMatrix(faceNormal);
    Eigen::RowVector2d projPoint = intersection * T;
    Eigen::MatrixX2d projPoly = (*polygonPtr) * T;

    if (isPointInPolygon2D(projPoly, projPoint, checkBoundary, false)) {
        intersections.emplace_back(intersection);
        return true;
    }
    else {
        return false;
    }
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if a 3D linear primitives (segments, rays, or lines) intersects with a 3D polygon.
 *
 * The function works with `Segment3D`, `Ray3D`, or `Line3D` types, all of which must have `P1` and `P2` as members of type
 * Eigen::RowVector3d representing the start and end points of the linear primitives.
 *
 * @tparam LineT: The type of the 3D linear primitive; it must provide `Eigen::RowVector3d` members `P1` and `P2`.
 * @param[in] polygonVertices: A matrix of size `N×3` where each row represents a 3D vertex of the 3D polygon.
 *                             The vertices may not be ordered or closed. The function will automatically close and sort the polygon if necessary.
 * @param[in] geoAB: The 3D linear primitive to check for intersection with the polygon.
 * @param[in] checkBoundary: If `true` (default), the function considers intersections on the polygon's boundaries.
 *                            If `false`, only intersections strictly inside the polygon are considered.
 * @param[in] needClosePolygon: If `true` (default), the function will automatically close the polygon by sorting its vertices
 *                             counter-clockwise and ensuring the polygon is closed.
 *                             If `false`, you must ensure that the `polygonVertices` is sorted and closed.
 *                             Setting this to `false` avoids the overhead of sorting and closing, improving performance when the polygon is already sorted and closed.
 * @return Returns `true` if at least one intersection is found; `false` otherwise.
 */
template<typename LineT>
bool isLinePolygonIntersection3D(
    const Eigen::MatrixX3d& polygonVertices,
    const LineT& geoAB,
    bool checkBoundary = true,
    bool needClosePolygon = true
) {
    if (!isCoplanar(polygonVertices)) return false;

    Eigen::MatrixX3d polygon;
    const Eigen::MatrixX3d* polygonPtr;
    if (needClosePolygon) {
        polygon = generateClosedPolygon(polygonVertices);
        polygonPtr = &polygon;
    }
    else {
        polygonPtr = &polygonVertices;
    }

    const Eigen::RowVector3d& refPoint = polygonPtr->row(0);
    Eigen::RowVector3d faceNormal = computePlaneNormal(*polygonPtr).first;

    const Eigen::RowVector3d dirAB = geoAB.P2 - geoAB.P1;
    const double denom = faceNormal.dot(dirAB);
    const Eigen::RowVector3d deltaA = geoAB.P1 - refPoint;
    const bool isParallel = std::abs(denom) < epsilon;
    if (checkBoundary && isParallel) {
        const Eigen::RowVector3d deltaB = geoAB.P2 - refPoint;
        if (std::abs(faceNormal.dot(deltaA)) > epsilon ||
            std::abs(faceNormal.dot(deltaB)) > epsilon) return false;

        Eigen::Matrix<double, 3, 2> T = compute3Dto2DTransformMatrix(faceNormal);
        Eigen::MatrixX3d deltaVertices = polygonPtr->rowwise() - refPoint;

        using Line2DType = typename line_traits<LineT>::Line2DType;
        if (isLinePolygonIntersection2D(deltaVertices * T, Line2DType{ deltaA * T, deltaB * T }, false)) {
            return true;
        }
        else {
            return false;
        }
    }

    if (isParallel) return false;

    const double t = faceNormal.dot(-deltaA) / denom;
    if (!line_traits<LineT>::validT(t)) return false;

    Eigen::RowVector3d intersection = geoAB.P1 + t * dirAB;
    Eigen::Matrix<double, 3, 2> T = compute3Dto2DTransformMatrix(faceNormal);
    Eigen::RowVector2d projPoint = intersection * T;
    Eigen::MatrixX2d projPoly = (*polygonPtr) * T;

    if (isPointInPolygon2D(projPoly, projPoint, checkBoundary, false)) {
        return true;
    }
    else {
        return false;
    }
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if a 3D linear primitives (segments, rays, or lines) intersects with a polyhedron and calculates intersection points.
 *
 * The function works with `Segment3D`, `Ray3D`, or `Line3D` types, all of which must have `P1` and `P2` as members of type
 * Eigen::RowVector3d representing the start and end points of the linear primitives.
 *
 * @tparam LineT: The type of the 3D linear primitive; it must provide `Eigen::RowVector3d` members `P1` and `P2`.
 * @param[in] polyhedronVertices: A matrix of size `N×3` where each row represents a 3D vertex of the polyhedron.
 * @param[in] geoAB: The 3D linear primitive to check for intersection with the polygon.
 * @param[out] intersections: If the lines intersect the polyhedron, this parameter will hold the intersection points, each of type `Eigen::RowVector3d`.
 * @return Returns `true` if at least one intersection is found; `false` otherwise.
 */
template<typename LineT>
bool isLinePolyhedronIntersection(
    const Eigen::MatrixX3d& polyhedronVertices,
    const LineT& geoAB,
    std::vector<Eigen::RowVector3d>& intersections
) {
    const Eigen::RowVector3d dirAB = geoAB.P2 - geoAB.P1;
    Mesh poly = buildConvexHullMesh(polyhedronVertices);
    for (auto face : poly.faces()) {
        auto he = poly.halfedge(face);
        auto start = he;
        Eigen::Matrix<double, 4, 3> eigenPoly;
        Eigen::Index i = 0;
        do {
            const Point_3 point = poly.point(poly.target(he));
            eigenPoly.row(i) << point.x(), point.y(), point.z();
            he = poly.next(he);
            i++;
        } while (he != start);
        eigenPoly.row(3) = eigenPoly.row(0);

        Eigen::RowVector3d faceNormal = computePlaneNormal(eigenPoly).first;
        double denom = faceNormal.dot(dirAB);
        if (std::abs(denom) < epsilon) continue;

        //const double D = faceNormal.dot(eigenPoly.row(0));
        const double t = (faceNormal.dot(eigenPoly.row(0)) - faceNormal.dot(geoAB.P1)) / denom;

        if (!line_traits<LineT>::validT(t)) continue;

        Eigen::RowVector3d intersection = geoAB.P1 + t * dirAB;
        Eigen::Matrix<double, 3, 2> T = compute3Dto2DTransformMatrix(faceNormal);

        if (isPointInPolygon2D(eigenPoly * T, intersection * T, false, false)) {
            intersections.emplace_back(intersection);
        }
    }
    return !intersections.empty();
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if a 3D linear primitives (segments, rays, or lines) intersects with a polyhedron and calculates intersection points.
 *
 * The function works with `Segment3D`, `Ray3D`, or `Line3D` types, all of which must have `P1` and `P2` as members of type
 * Eigen::RowVector3d representing the start and end points of the linear primitives.
 *
 * @tparam LineT: The type of the 3D linear primitive; it must provide `Eigen::RowVector3d` members `P1` and `P2`.
 * @param[in] polyhedronVertices: A matrix of size `N×3` where each row represents a 3D vertex of the polyhedron.
 * @param[in] geoAB: The 3D linear primitive to check for intersection with the polygon.
 * @return Returns `true` if at least one intersection is found; `false` otherwise.
 */
template<typename LineT>
bool isLinePolyhedronIntersection(
    const Eigen::MatrixX3d& polyhedronVertices,
    const LineT& geoAB
) {
    const Eigen::RowVector3d dirAB = geoAB.P2 - geoAB.P1;
    Mesh poly = buildConvexHullMesh(polyhedronVertices);
    for (auto face : poly.faces()) {
        auto he = poly.halfedge(face);
        auto start = he;
        Eigen::Matrix<double, 4, 3> eigenPoly;
        Eigen::Index i = 0;
        do {
            const Point_3 point = poly.point(poly.target(he));
            eigenPoly.row(i) << point.x(), point.y(), point.z();
            he = poly.next(he);
            i++;
        } while (he != start);
        eigenPoly.row(3) = eigenPoly.row(0);

        Eigen::RowVector3d faceNormal = computePlaneNormal(eigenPoly).first;
        double denom = faceNormal.dot(dirAB); 
        if (std::abs(denom) < epsilon) continue;

        //const double D = faceNormal.dot(eigenPoly.row(0));
        const double t = (faceNormal.dot(eigenPoly.row(0)) - faceNormal.dot(geoAB.P1)) / denom;

        if (!line_traits<LineT>::validT(t)) continue;

        Eigen::RowVector3d intersection = geoAB.P1 + t * dirAB;
        Eigen::Matrix<double, 3, 2> T = compute3Dto2DTransformMatrix(faceNormal);

        if (isPointInPolygon2D(eigenPoly * T, intersection * T, false, false)) return true;

    }
    return false;
}
//************************************************************************************************************************//

#endif