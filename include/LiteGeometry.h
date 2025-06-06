#ifndef LITEGEOMETRY_H_
#define LITEGEOMETRY_H_
#include "Utility.h"
#include "PointTool.h"
#include "LineTool.h"
////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Filters a set of 2D points based on their location relative to a polygon.
 *
 * This function filters a set of points by checking if each point lies inside, outside, or on the boundary of a polygon.
 * The behavior of the filter is controlled by the `removeBoundary` and `removeInside` flags:
 * - When `removeInside` is false and `removeBoundary` is false, only external points are removed, keeping internal and boundary points.
 * - When `removeInside` is false and `removeBoundary` is true, both external points and boundary points are removed, keeping only internal points.
 * - When `removeInside` is true and `removeBoundary` is false, only internal points are removed, keeping external and boundary points.
 * - When `removeInside` is true and `removeBoundary` is true, both internal and boundary points are removed, keeping only external points.
 *
 * @param[in] polygonVertices: A matrix (`Eigen::MatrixX2d`) representing the vertices of the polygon in 2D space.
 * @param[in] points: A vector of 2D points (`std::vector<Eigen::RowVector2d>`) to be filtered.
 * @param[in] removeBoundary: A boolean flag to control whether boundary points are removed (Default: true).
 * @param[in] removeInside: A boolean flag to control whether inside points are removed (Default: true).
 * @param[in] needClosePolygon: A boolean flag to control whether the polygon vertices should be used to generate a closed polygon before filtering (Default: true).
 * @return A vector (`std::vector<Eigen::RowVector2d>`) containing the filtered points.
 */
std::vector<Eigen::RowVector2d> filterPointsByPolygon(
    const Eigen::MatrixX2d& polygonVertices,
    const std::vector<Eigen::RowVector2d>& points,
    bool removeBoundary = true,
    bool removeInside = true,
    bool needClosePolygon = true
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Filters a set of 3D points based on their location relative to a polyhedron.
 *
 * This function filters a set of 3D points by checking whether each point lies inside, outside, or on the boundary of a polyhedron.
 * The behavior of the filter is controlled by the `removeBoundary` and `removeInside` flags:
 * - When `removeInside` is false and `removeBoundary` is false, only external points are removed, keeping internal and boundary points.
 * - When `removeInside` is false and `removeBoundary` is true, both external points and boundary points are removed, keeping only internal points.
 * - When `removeInside` is true and `removeBoundary` is false, only internal points are removed, keeping external and boundary points.
 * - When `removeInside` is true and `removeBoundary` is true, both internal and boundary points are removed, keeping only external points.
 *
 * @param[in] polyhedronVertices: A matrix (`Eigen::MatrixX3d`) representing the vertices of the polyhedron in 3D space.
 * @param[in] points: A vector of 3D points (`std::vector<Eigen::RowVector2d>`) to be filtered.
 * @param[in] removeBoundary: A boolean flag to control whether boundary points are removed (Default: true).
 * @param[in] removeInside: A boolean flag to control whether inside points are removed (Default: true).
 * @return A vector (`std::vector<Eigen::RowVector3d>`) containing the filtered points.
 */
std::vector<Eigen::RowVector3d> filterPointsByPolyhedron(
    const Eigen::MatrixX3d& polyhedronVertices,
    const std::vector<Eigen::RowVector3d>& points,
    bool removeBoundary = true,
    bool removeInside = true
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Computes a PCA-aligned rectangle for a set of 2D points.
 *
 * Calculates the rectangle oriented along the principal axes of the input points.
 *
 * @param[in] points: A matrix of 2D points (`Eigen::MatrixX2d`, each row is a point).
 * @param[in] useConvexHull: A boolean flag indicating whether to compute the convex hull of the points (Default: false).
 * @return A pair consisting of:
 *         - A 4x2 matrix (`Eigen::MatrixX2d`) representing the corners of the PCA rectangle.
 *         - A 2D vector representing the centroid (mean) of the rectangle.
 */
std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computePCAOBB2D(
    const Eigen::MatrixX2d& points,
    bool useConvexHull = false
);

/**
 * @brief Computes a PCA-aligned rectangle for a set of 2D points.
 *
 * Calculates the rectangle oriented along the principal axes of the input points.
 *
 * @param[in] points: A vector of 2D points (`std::vector<Eigen::RowVector2d>`, each representing a 2D point).
 * @param[in] useConvexHull: If true, computes the convex hull before PCA processing. (Default: false).
 * @return A pair consisting of:
 *         - A 4x2 matrix (`Eigen::Matrix<double, 4, 2>`) representing the corners of the PCA rectangle.
 *         - A 2D vector representing the centroid (mean) of the rectangle.
 */
std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computePCAOBB2D(
    const std::vector<Eigen::RowVector2d>& points,
    bool useConvexHull = false
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Computes a PCA-aligned box for a set of 3D points.
 *
 * Calculates the bounding box oriented along the principal axes of the input points.
 *
 * @param[in] points: A matrix of 3D points (`Eigen::MatrixX3d`, each row is a point).
 * @param[in] useConvexHull: A boolean flag indicating whether to compute the convex hull of the points.
 * @return A pair consisting of:
 *         - A 8x3 matrix (`Eigen::Matrix<double, 8, 3>`, each row is a vertex) representing the corners of the PCA box.
 *         - A 3D vector (`Eigen::RowVector2d`) representing the centroid of the box.
 */
std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computePCAOBB3D(
    const Eigen::MatrixX3d& points,
    bool useConvexHull = false
);

/**
 * @brief Computes a PCA-aligned box for a set of 3D points.
 *
 * Calculates the bounding box oriented along the principal axes of the input points.
 *
 * @param[in] points: A vector of 3D points (`std::vector<Eigen::RowVector3d>`).
 * @param[in] useConvexHull: A boolean flag indicating whether to compute the convex hull of the points.
 * @return A pair consisting of:
 *         - A 8x3 matrix (`Eigen::Matrix<double, 8, 3>`, each row is a vertex) representing the corners of the PCA box.
 *         - A 3D vector (`Eigen::RowVector3d`) representing the centroid of the box.
 */
std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computePCAOBB3D(
    const std::vector<Eigen::RowVector3d>& points,
    bool useConvexHull = false
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Computes the minimum bounding rectangle for a set of 2D points.
 *
 * Calculates the smallest oriented rectangle that encloses all the 2D points by applying a rotating calipers approach
 * on the convex hull. This method generates a tighter and more optimal bounding rectangle compared to PCA-based methods.
 * 
 * @param[in] points A matrix of 2D points (`Eigen::MatrixX2d`, each row is a point).
 * @return A pair consisting of:
 *         - A 4x2 matrix (`Eigen::Matrix<double, 4, 2>`, each row is a vertex) representing the corners of the minimum bounding rectangle.
 *         - A 2D vector (`Eigen::RowVector2d`) representing the centroid of the rectangle.
 */
std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computeMinOBB2D(
    const Eigen::MatrixX2d& points
);

/**
 * @brief Computes the minimum bounding rectangle for a set of 2D points.
 *
 * Calculates the smallest oriented rectangle that encloses all the 2D points by applying a rotating calipers approach
 * on the convex hull. This method generates a tighter and more optimal bounding rectangle compared to PCA-based methods.
 *
 * @param[in] points: A vector of 2D points (`std::vector<Eigen::RowVector2d>`).
 * @return A pair consisting of:
 *         - A 4x2 matrix (`Eigen::Matrix<double, 4, 2>`, each row is a vertex) representing the corners of the minimum bounding rectangle.
 *         - A 2D vector (`Eigen::RowVector2d`) representing the centroid of the rectangle.
 */
std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computeMinOBB2D(
    const std::vector<Eigen::RowVector2d>& points
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Computes the minimum bounding box for a set of 3D points.
 *
 * Calculates the smallest oriented box that encloses all the 3D points by applying geometric optimization
 * on the convex hull. This method generates a tighter and more optimal bounding box compared to PCA-based methods.
 * 
 * @param[in] points: A matrix of 3D points (`Eigen::MatrixX3d`, each row is a point).
 * @return A pair consisting of:
 *         - A 8x3 matrix (`Eigen::Matrix<double, 8, 3>`, each row is a vertex) representing the corners of the minimum bounding box.
 *         - A 3D vector (`Eigen::RowVector3d`) representing the centroid of the box.
 */
std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computeMinOBB3D(
    const Eigen::MatrixX3d& points
);

/**
 * @brief Computes the minimum bounding box for a set of 3D points.
 *
 * Calculates the smallest oriented box that encloses all the 3D points by applying geometric optimization
 * on the convex hull. This method generates a tighter and more optimal bounding box compared to PCA-based methods.
 *
 * @param[in] points: A vector of 3D points (`std::vector<Eigen::RowVector3d>`).
 * @return A pair consisting of:
 *         - A 8x3 matrix (`Eigen::Matrix<double, 8, 3>`, each row is a vertex) representing the corners of the minimum bounding box.
 *         - A 3D vector (`Eigen::RowVector3d`) representing the centroid of the box.
 */
std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computeMinOBB3D(
    const std::vector<Eigen::RowVector3d>& points
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Converts 3D world coordinates to 2D camera image coordinates.
 *
 * This function transforms a set of 3D target vertices from world coordinates to 2D camera image coordinates.
 * The transformation uses the camera position, the direction of the camera (from the camera position to the camera direction point),
 * and the focal length (FL) of the camera to convert the target vertices from world coordinates to 2D camera image coordinates.
 *
 * @param[in] targetVertices: A matrix of 3D target vertices (`Eigen::MatrixX3d`, each row represents a vertex).
 * @param[in] cameraDirPoint: A 3D point (`Eigen::RowVector3d`) representing the direction of the camera.
 * @param[in] cameraPosition: A 3D point (`Eigen::RowVector3d`) representing the position of the camera.
 * @param[in] FL: The focal length of the camera.
 * @return A matrix of 2D image coordinates (`Eigen::MatrixX2d`) corresponding to the input 3D target vertices.
 */
Eigen::MatrixX2d WorldToCameraImageCoords(
    const Eigen::MatrixX3d& targetVertices,
    const Eigen::RowVector3d& cameraDirPoint,
    const Eigen::RowVector3d& cameraPosition,
    double FL
);

/**
 * @brief Converts 3D world coordinates to 2D camera image coordinates (using the centroid as the camera direction point).
 *
 * This function is a simplified version of the previous function where the camera direction point is assumed to be
 * the centroid (mean) of the target vertices. It uses the camera position and focal length (FL) to convert the target
 * vertices from world coordinates to 2D camera image coordinates.
 *
 * @param[in] targetVertices A matrix of 3D target vertices (`Eigen::MatrixX3d`, each row represents a vertex).
 * @param[in] cameraPosition: A 3D point (`Eigen::RowVector3d`) representing the position of the camera.
 * @param[in] FL: The focal length of the camera.
 * @return A matrix of 2D image coordinates (`Eigen::MatrixX2d`) corresponding to the input 3D target vertices.
 */
Eigen::MatrixX2d WorldToCameraImageCoords(
    const Eigen::MatrixX3d& targetVertices,
    const Eigen::RowVector3d& cameraPosition,
    double FL
);
//************************************************************************************************************************//

#endif