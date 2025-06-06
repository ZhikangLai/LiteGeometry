#ifndef RASTERIZAYIONTOOL_H
#define RASTERIZAYIONTOOL_H
#include "AbslEigenRowVectoriHash.h"

/**
 * @brief Rasterizes a 2D line segment and its neighboring cells using Bresenham's algorithm.
 *
 * This function applies Bresenham's line algorithm to compute discrete 2D grid points
 * that approximate a line between two points `p1` and `p2`.
 *
 * The points are stored as a unique set of 2D integer coordinates in `RowVector2iSet`,
 * which is implemented using `absl::flat_hash_set` for optimal performance.
 * 
 * @param[in] p1: A 2D point representing the start of the line.
 * @param[in] p2: A 2D point representing the end of the line.
 * @return A set of 2D integer grid points including the main line and its offset neighbors.
 */
RowVector2iSet rasterizeBresenhamLine2D(
    const Eigen::RowVector2d& p1,
    const Eigen::RowVector2d& p2
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Rasterizes a 2D line segment and its neighboring cells using Bresenham's algorithm.
 *
 * This function applies Bresenham's line algorithm to compute discrete 2D grid points
 * that approximate a line between two points `p1` and `p2`.
 *
 * The points are stored as a unique set of 2D integer coordinates in `RowVector2iSet`,
 * which is implemented using `absl::flat_hash_set` for optimal performance.
 * 
 * @param[in] p1: A 2D point representing the start of the line.
 * @param[in] p2: A 2D point representing the end of the line.
 * @param[out] rawLinePoints: A vector of `Eigen::RowVector2i` storing the raw rasterized line cells before offset expansion.
 * @return A set of 2D integer grid points including the main line and its offset neighbors.
 */
RowVector2iSet rasterizeBresenhamLine2D(
    const Eigen::RowVector2d& p1,
    const Eigen::RowVector2d& p2,
    std::vector<Eigen::RowVector2i>& rawLinePoints
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Rasterizes a 2D polygon (including its interior and surrounding grid cells).
 *
 * This function rasterizes a 2D polygon by:
 * 1. Rasterizing its boundary edges using Bresenham's line algorithm.
 * 2. Filling the interior grid points by testing whether each grid cell center lies inside the polygon.
 * 3. Expanding the result by including offset neighbors for each grid cell to ensure a wider coverage.
 *
 * The points are stored as a unique set of 2D integer coordinates in `RowVector2iSet`,
 * which is implemented using `absl::flat_hash_set` for optimal performance.
 *
 * @param[in] polygonVertices: A matrix of 2D points representing the vertices of the polygon.
 * @param[in] needClosePolygon: If `true` (default), the function will automatically close the 2D polygon by sorting its vertices
 *                              counter-clockwise and ensuring the polygon is closed.
 *                              If `false`, you must ensure that the `polygonVertices` is sorted and closed.
 *                              Setting this to `false` avoids the overhead of sorting and closing, improving performance when the polygon is already sorted and closed.
 * @return A set of 2D grid points covering the polygon's boundary, interior, and surrounding cells.
 */
RowVector2iSet rasterizePolygon2D(
    const Eigen::MatrixX2d& polygonVertices,
    bool needClosePolygon = true
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Rasterizes a 3D line segment and its neighboring cells using Bresenham's algorithm.
 *
 * This function applies Bresenham's line algorithm to compute discrete 3D grid points
 * that approximate a line between two points `p1` and `p2`.
 *
 * The grid points are stored in a unique set of 3D integer coordinates in  `RowVector3iSet`, defined as:
 * using RowVector3iSet = absl::flat_hash_set<Eigen::RowVector3i, EigenRowVectoriHash, EigenRowVectoriEqual>;
 *
 * @param[in] p1: A 3D point (Eigen::RowVector3d) representing the start of the line.
 * @param[in] p2: A 3D point (Eigen::RowVector3d) representing the end of the line.
 * @return A unique set (RowVector3iSet) of 3D grid cells that includes both the raw cells and their neighboring cells.
 */
RowVector3iSet rasterizeBresenhamLine3D(
    const Eigen::RowVector3d& p1,
    const Eigen::RowVector3d& p2
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Rasterizes a 3D line segment and its neighboring cells using Bresenham's algorithm.
 * 
 * This function applies Bresenham's algorithm to compute discrete 3D grid points that best approximate
 * a line segment from `p1` to `p2`. In addition, the function outputs the raw grid cells 
 * (i.e. the primary cells before neighboring expansion) into the provided vector `rawLinePoints`. 
 *
 * The grid points are stored as a unique set of 3D integer coordinates in  `RowVector3iSet`, defined as:
 * using RowVector3iSet = absl::flat_hash_set<Eigen::RowVector3i, EigenRowVectoriHash, EigenRowVectoriEqual>;
 * 
 * @param[in] p1: A 3D point (Eigen::RowVector3d) representing the start of the line.
 * @param[in] p2: A 3D point (Eigen::RowVector3d) representing the end of the line.
 * @param[out] rawLinePoints: If provided, this parameter will hold the raw grid cells 
 *             (std::vector<Eigen::RowVector3i>) computed by Bresenham's algorithm before neighboring expansion.
 * @return A unique set (RowVector3iSet) of 3D grid cells that includes both the raw cells and their neighboring cells.
 * 
 */
RowVector3iSet rasterizeBresenhamLine3D(
    const Eigen::RowVector3d& p1,
    const Eigen::RowVector3d& p2,
    std::vector<Eigen::RowVector3i>& rawLinePoints
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Rasterizes a 3D polygon (including its interior and surrounding grid cells).
 *
 * This function rasterizes a 3D polygon by:
 * 1. Rasterizing its boundary edges using Bresenham's line algorithm.
 * 2. Filling the interior grid points by testing whether each grid cell center lies inside the polygon.
 * 3. Expanding the result by including offset neighbors for each grid cell to ensure a wider coverage.
 * 
 * The grid points are stored as a unique set of 3D integer coordinates in a `RowVector3iSet`, defined as:
 * using RowVector3iSet = absl::flat_hash_set<Eigen::RowVector3i, EigenRowVectoriHash, EigenRowVectoriEqual>.
 * 
 * 
 * @param[in] polygonVertices: A matrix of 3D points (Eigen::MatrixX3d) representing the vertices of the polygon.
 * @param[in] needClosePolygon: If `true` (default), the function will automatically close the 3D polygon by sorting its vertices
 *                              counter-clockwise and ensuring the polygon is closed.
 *                              If `false`, you must ensure that the `polygonVertices` is sorted and closed.
 *                              Setting this to `false` avoids the overhead of sorting and closing, improving performance when the polygon is already sorted and closed.
 * @return A set of 3D grid points covering the polygon's boundary, interior, and surrounding cells.
 */
RowVector3iSet rasterizePolygon3D(
    const Eigen::MatrixX3d& polygonVertices,
    bool needClosePolygon = true
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Rasterizes a 3D polyhedron (including its faces and surrounding grid cells).
 *
 * This function rasterizes a 3D polyhedron by:
 * 1. Rasterizing its boundary edges using the `rasterizeBresenhamLine3D` algorithm.
 * 2. Generating grid points on the faces of the polyhedron by rasterizing each polygonal face's boundary using the `rasterizePolygon3D` algorithm.
 * 3. Expanding the result by including offset neighbors for each grid cell to ensure a wider coverage.
 *
 * The grid points are stored as a unique set of 3D integer coordinates in a `RowVector3iSet`, defined as:
 * using RowVector3iSet = absl::flat_hash_set<Eigen::RowVector3i, EigenRowVectoriHash, EigenRowVectoriEqual>.
 *
 * @param[in] polyhedronVertices: A matrix of 3D points (Eigen::MatrixX3d) representing the vertices of the polyhedron.

 * @return A unique set (RowVector3iSet) of 3D grid points representing the faces of the polyhedron and the surrounding cells.
 */
RowVector3iSet rasterizePolyhedron(
    const Eigen::MatrixX3d& polyhedronVertices
);
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Rasterizes a 3D polyhedron (including its faces and surrounding grid cells).
 *
 * This function rasterizes a 3D polyhedron by:
 * 1. Rasterizing its boundary edges using the `rasterizeBresenhamLine3D` algorithm.
 * 2. Generating grid points on the faces of the polyhedron by rasterizing each polygonal face's boundary using the `rasterizePolygon3D` algorithm.
 * 3. Expanding the result by including offset neighbors for each grid cell to ensure a wider coverage.
 *
 * The grid points are stored as a unique set of 3D integer coordinates in a `RowVector3iSet`, defined as:
 * using RowVector3iSet = absl::flat_hash_set<Eigen::RowVector3i, EigenRowVectoriHash, EigenRowVectoriEqual>.
 *
 * @param[in] polyhedronMesh: A Mesh object representing the polyhedron. The mesh is expected to define all faces of the polyhedron,
 *                           from which vertices for each face are extracted for rasterization.
 *
 * @return A unique set (RowVector3iSet) of 3D grid points representing the faces of the polyhedron and the surrounding cells.
 */
RowVector3iSet rasterizePolyhedron(
    const Mesh& polyhedronMesh
);
//************************************************************************************************************************//
#endif
