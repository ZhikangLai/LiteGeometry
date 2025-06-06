#ifndef LITEGEOMETRYBASE_H_
#define LITEGEOMETRYBASE_H_
#include "LineTraits.h"

/**
 * @brief Builds a 3D convex hull mesh from a set of polyhedron vertices.

 * @param[in] polyhedronVertices: A matrix of size N¡Á3 representing the vertices of the polyhedron.
 *                               Each row represents a 3D point (x, y, z).
 * @return A `Mesh` object representing the convex hull of the input vertices..
 */
Mesh buildConvexHullMesh(
    const Eigen::MatrixX3d& polyhedronVertices
);

/**
 * @brief Generates a closed polygon from a set of unordered 2D or coplanar 3D points.
 *
 * This function computes the centroid of the input points, sorts them counter-clockwise
 * based on the angle relative to the centroid, and appends the first point to the end
 * if the polygon is not already closed.
 *
 * Supports both 2D (Nx2) and 3D (Nx3) point inputs. For 3D points, coplanarity is required.
 *
 * @param[in] polygonVertices: A matrix of size N¡Á2 or N¡Á3 representing the polygon's vertices.
 *
 * @return A matrix (`Eigen::MatrixXd`) of sorted polygon vertices in counter-clockwise order,
 *         closed by duplicating the first vertex at the end if necessary.
 *
 */
Eigen::MatrixXd generateClosedPolygon(
    const Eigen::MatrixXd& polygonVertices
);

/**
 * @brief Computes the normal vector(s) of a 3D plane defined by a set of coplanar points.
 *
 * Given at least three 3D points that lie on the same plane, this function computes the plane's unit normal vector
 * using the cross product of two edges formed by the first three points. Since a plane has two opposite normal directions,
 * both are returned.
 *
 * @param[in] planePoints: A matrix of size `N¡Á3` where each row represents a 3D point on the same plane.
 *                         At least 3 points are required, and all points must be coplanar.
 *
 * @return A `std::pair<Eigen::RowVector3d, Eigen::RowVector3d>` containing the two opposite unit normal vectors of the plane.
 *         The first is the computed normal, the second is its negation (i.e., the opposite direction).
 *
 */
std::pair<Eigen::RowVector3d, Eigen::RowVector3d> computePlaneNormal(
    const Eigen::MatrixX3d& planePoints
);


/*
 * @brief Computes a 3D-to-2D transformation matrix based on the given plane normal.
 *
 * This function generates a 3¡Á2 matrix whose columns form an orthonormal basis lying on the
 * plane that is perpendicular to the specified normal vector. The resulting matrix can be used
 * to map 3D points into 2D coordinates on this plane.
 *
 * @param[in] planeNormal: The normal vector of the target plane (Eigen::RowVector3d).
 *
 * @return Return a `Eigen::Matrix<double, 3, 2>` matrix whose columns are orthonormal vectors
 *         spanning the plane orthogonal to the input normal.
 *
 */
Eigen::Matrix<double, 3, 2> compute3Dto2DTransformMatrix(
    const Eigen::RowVector3d& planeNormal
);


/**
 * @brief Checks if a set of 3D points are coplanar.
 *
 * This function checks whether the provided set of 3D points lies on the same plane.
 * At least 3 points (i.e., a matrix with at least 3 rows) are required for a valid coplanarity check.
 *
 * @param[in] points: A matrix of size N¡Á3, where each row represents a 3D point (x, y, z).
 * @return `true` if the points are coplanar, `false` otherwise.
 */
bool isCoplanar(const Eigen::MatrixX3d& points);

#endif