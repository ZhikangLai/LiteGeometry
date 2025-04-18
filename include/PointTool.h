#ifndef POINTTOOL_H_
#define POINTTOOL_H_
#include "LineTraits.h"
#include "LightGeoBase.h"

//*****************************************************************************************//
/**
 * @brief Checks if a 2D point lies on a given 2D linear primitive.
 *
 * The function works with `Segment2D`, `Ray2D`, or `Line2D` types, all of which must have `P1` and `P2` as members of type
 * `Eigen::RowVector2d` representing the start and end points of the linear primitives.
 *
 * @tparam LineT: The type of the 2D linear primitive; it must provide `Eigen::RowVector2d` members `P1` and `P2`.
 * @param[in] geoAB: The geometric primitive to test.
 * @param[in] point: The 2D point (Eigen::RowVector2d) that is being tested.
 * 
 * @return Returns `true` if the point lies on the linear primitive; `false` otherwise.
 */
template<typename LineT>
bool isPointOnLine2D(
    const LineT& geoAB,
    const Eigen::RowVector2d& point
) {
    const Eigen::RowVector2d dirAB = geoAB.P2 - geoAB.P1;
    const Eigen::RowVector2d dirAP = point - geoAB.P1;
    const double lenSqAB = dirAB.squaredNorm();

    if (lenSqAB < highEpsilon) return dirAP.squaredNorm() < highEpsilon;

    const double crossArea = dirAB.cross(dirAP);
    const double crossSquaredNorm = crossArea * crossArea;
    if (crossSquaredNorm > highEpsilon * lenSqAB) return false;

    const double tScaled = dirAP.dot(dirAB);
    return line_traits<LineT>::validTScaled(tScaled, lenSqAB);
};
//*****************************************************************************************//

//*****************************************************************************************//
/**
 * @brief Checks if a 3D point lies on a given 3D linear primitive.
 *
 * The function works with `Segment3D`, `Ray3D`, or `Line3D` types, all of which must have `P1` and `P2` as members of type
 * `Eigen::RowVector3d` representing the start and end points of the linear primitives.
 *
 * @tparam LineT: The type of the 3D linear primitive; it must provide `Eigen::RowVector3d` members `P1` and `P2`.
 * @param[in] geoAB: The geometric primitive to test.
 * @param[in] point: The 3D point (Eigen::RowVector3d) that is being tested.
 * 
 * @return Returns `true` if the point lies on the linear primitive; `false` otherwise.
 */
template<typename LineT>
bool isPointOnLine3D(
    const LineT& geoAB,
    const Eigen::RowVector3d& point
) {

    const Eigen::RowVector3d dirAB = geoAB.P2 - geoAB.P1;
    const Eigen::RowVector3d dirAP = point - geoAB.P1;

    const double lenSqAB = dirAB.squaredNorm();
    if (lenSqAB < highEpsilon) {
        return (dirAP.squaredNorm() < highEpsilon);
    }

    const double crossSquaredNorm = dirAB.cross(dirAP).squaredNorm();

    if (crossSquaredNorm > highEpsilon * lenSqAB) return false;

    const double tScaled = dirAP.dot(dirAB);

    return line_traits<LineT>::validTScaled(tScaled, lenSqAB);
};
//*****************************************************************************************//

//*****************************************************************************************//
/**
 * @brief Checks if a 2D point lies inside a 2D polygon or on its boundary.
 *
 * @param[in] polygonVertices: A matrix of size `N¡Á2` where each row represents a 2D vertex of the polygon.
 * @param[in] point: The 3D point (`Eigen::RowVector2d`) that is being tested.
 * @param[in] checkBoundary: If `true` (default), the function considers intersections on the polygon's boundaries.
 *                            If `false`, only intersections strictly inside the polygon are considered.
 * @param[in] needClosePolygon: If `true` (default), the function will automatically close the 2D polygon by sorting its vertices
 *                              counter-clockwise and ensuring the polygon is closed.
 *                              If `false`, you must ensure that the `polygonVertices` is sorted and closed.
 *                              Setting this to `false` avoids the overhead of sorting and closing, improving performance when the polygon is already sorted and closed.
 *
 * @return `true` if the point lies inside the polygon or on the boundary (depending on `checkBoundary`), `false` otherwise.
 *
 */
bool isPointInPolygon2D(
    const Eigen::MatrixX2d& polygonVertices,
    const Eigen::RowVector2d& point,
    bool checkBoundary = true,
    bool needClosePolygon = true
);
//*****************************************************************************************//

//*****************************************************************************************//
/**
 * @brief Checks if a 3D point lies inside a 3D polygon or on its boundary.
 *
 * @param[in] polygonVertices: A matrix of size `N¡Á3` where each row represents a 3D vertex of the polygon.
 * @param[in] point: The 3D point (`Eigen::RowVector3d`) that is being tested.
 * @param[in] checkBoundary: If `true` (default), the function considers intersections on the polygon's boundaries.
 *                           If `false`, only intersections strictly inside the polygon are considered.
 * @param[in] needClosePolygon: If `true` (default), the function will automatically close the 3D polygon by sorting its vertices
 *                              counter-clockwise and ensuring the polygon is closed.
 *                              If `false`, you must ensure that the `polygonVertices` is sorted and closed.
 *                              Setting this to `false` avoids the overhead of sorting and closing, improving performance when the polygon is already sorted and closed.
 *
 * @return `true` if the point lies inside the polygon or on the boundary (depending on `checkBoundary`), `false` otherwise.
 *
 */
bool isPointInPolygon3D(
    const Eigen::MatrixX3d& polygonVertices,
    const Eigen::RowVector3d& point,
    bool checkBoundary = true,
    bool needClosePolygon = true
);
//*****************************************************************************************//

//*****************************************************************************************//
/**
 * @brief Checks if a 3D point lies inside a polyhedron or on its boundary.
 *
 * @param[in] polyhedronVertices: A matrix of size `N¡Á3` where each row represents a 3D vertex of the polyhedron.
 * @param[in] point: The 3D point (`Eigen::RowVector3d`) that is being tested.
 * @param[in] checkBoundary: If `true` (default), the function considers intersections on the polygon's boundaries.
 *                           If `false`, only intersections strictly inside the polygon are considered.
 * 
 * @return `true` if the point lies inside the polyhedron or on the boundary (depending on `checkBoundary`), `false` otherwise.
 */
bool isPointInPolyhedron(
    const Eigen::MatrixX3d& polyhedronVertices,
    const Eigen::RowVector3d& point,
    bool checkBoundary = true
);
//*****************************************************************************************//
#endif