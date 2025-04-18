#include "lightGeo.h"
#include <execution>
#include <numeric> 


/**
 * @brief Builds a 3D convex hull mesh from a set of polyhedron vertices.
 *
 * This function takes a set of 3D vertices as input, computes the convex hull using CGAL's
 * `convex_hull_3` algorithm, and returns the resulting convex hull mesh.
 * The input points must represent vertices of a polyhedron in 3D space.
 *
 * @param[in] polyhedronVertices A matrix of size N×3 representing the vertices of the polyhedron.
 *                               Each row represents a 3D point (x, y, z).
 * @return Mesh The mesh representing the convex hull of the input vertices.
 */
Mesh buildConvexHullMesh(const Eigen::MatrixX3d& polyhedronVertices) {
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


static Eigen::MatrixXd computeCovarianceMatrix(
    const Eigen::MatrixXd& data,
    bool unbiased = true
) {
    Eigen::MatrixXd centered = data.rowwise() - data.colwise().mean();
    double divisor = unbiased ? static_cast<double>(data.rows() - 1) : static_cast<double>(data.rows());
    return (centered.transpose() * centered) / divisor;
}

//************************************************************************************************************************//
/**
 * @brief Generates a closed polygon from a set of unordered 2D or coplanar 3D points.
 *
 * This function computes the centroid of the input points, sorts them counter-clockwise
 * based on the angle relative to the centroid, and appends the first point to the end
 * if the polygon is not already closed.
 *
 * Supports both 2D (Nx2) and 3D (Nx3) point inputs. For 3D points, coplanarity is required.
 *
 * @param[in] polygonVertices A matrix of size N×2 or N×3 representing the polygon's vertices.
 *
 * @return Eigen::MatrixXd A matrix of polygon vertices arranged in counter-clockwise order,
 *                         with the first vertex duplicated at the end if the polygon is not already closed.
 *
 */
Eigen::MatrixXd generateClosedPolygon(const Eigen::MatrixXd& polygonVertices) {
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
/**
 * @brief Computes the normal vector(s) of a plane defined by a set of 3D coplanar points.
 *
 * This function verifies that at least three points are provided and that they are coplanar.
 * It then computes a unit normal vector of the plane using the cross product of two vectors
 * derived from the points. Both the normal and its opposite direction are returned.
 *
 * @param[in] planePoints A matrix of size n×3 (Eigen::MatrixX3d), each row representing a 3D point.
 *                        At least three coplanar points are required.
 * @return std::pair<Eigen::RowVector3d, Eigen::RowVector3d> A pair of unit normal vectors:
 *         the computed normal and its opposite direction.
 *
 * @throws std::invalid_argument if fewer than three points are provided or if the points are not coplanar.
 */
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
/*
 * @brief Computes a projection matrix that maps 3D vectors onto a 2D plane defined by its normal.
 *
 * This function generates a 3×2 matrix whose columns form an orthonormal basis
 * lying on the specified plane. The matrix can be used to project 3D points onto this 2D plane.
 *
 * @param[in] planeNormal The normal vector of the target plane (Eigen::RowVector3d, 1×3 row vector).
 *
 * @return Eigen::Matrix<double, 3, 2> A 3×2 matrix whose columns are orthonormal vectors
 *         spanning the plane orthogonal to the input normal.
 *
 */
Eigen::Matrix<double, 3, 2> computeProjectionMatrix(const Eigen::RowVector3d& planeNormal) {
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
/**
 * @brief Checks if a set of 3D points are coplanar.
 *
 * This function checks whether the provided set of 3D points lies on the same plane.
 * At least 3 points (i.e., a matrix with at least 3 rows) are required for a valid coplanarity check.
 *
 * @param[in] points A matrix of size N×3, where each row represents a 3D point (x, y, z).
 * @return `true` if the points are coplanar, `false` otherwise.
 */
bool isCoplanar(const Eigen::MatrixX3d& points) {
    Eigen::MatrixX3d centered = points.bottomRows(points.rows() - 1).rowwise() - points.row(0);
    Eigen::JacobiSVD<Eigen::MatrixX3d> svd(centered, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Eigen::VectorXd& singularValues = svd.singularValues();
    return singularValues.size() < 3 || singularValues(2) <= epsilon;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if a point lies on a 2D line.
 *
 * This function determines if a given point lies on a 2D line, based on the specified `lineType`.
 * The line can be interpreted as a segment, a ray, or an infinite line:
 * - If `lineType` is 1, it checks if the point lies on the **line segment** from A to B.
 * - If `lineType` is 2, it checks if the point lies on the **ray** starting at A and going through B.
 * - If `lineType` is 3, it checks if the point lies on the **infinite line** defined by A and B.
 *
 * @param[in] A The starting point of the line (1×2 vector).
 * @param[in] B The ending point of the line (1×2 vector).
 * @param[in] point The point to check (1×2 vector).
 * @param[in] lineType Specifies the type of line:
 *                     1 = Line segment,
 *                     2 = Ray extending from point A,
 *                     3 = Infinite line. (default: 1)
 * @return `true` if the point lies on the line based on `lineType`, `false` otherwise.
 */
bool isPointOnLine2D(
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    const Eigen::RowVector2d& point,
    int lineType
) {
    if (lineType < 1 || lineType > 3) {
        throw std::invalid_argument("Invalid lineType: must be 1 (segment), 2 (ray), or 3 (line)");
    }
    const Eigen::RowVector2d AB = B - A;
    const Eigen::RowVector2d AP = point - A;
    const double lengthSquaredAB = AB.squaredNorm();

    if (lengthSquaredAB < highEpsilon) return AP.squaredNorm() < highEpsilon;

    const double cross = AB.cross(AP);
    const double crossSquared = cross * cross;
    if (crossSquared > highEpsilon * lengthSquaredAB) return false;

    switch (lineType) {
    case 1: {
        const double AP_dot_AB = AP.dot(AB);
        return (AP_dot_AB >= -epsilon) && (AP_dot_AB <= lengthSquaredAB + epsilon);
    }
    case 2: {
        const double AP_dot_AB = AP.dot(AB);
        return AP_dot_AB >= -epsilon;
    }
    case 3:
        return true;
    }
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if a point lies on a 3D line.
 *
 * This function determines if a given point lies on a 3D line, based on the specified `lineType`.
 * The line can be interpreted as a segment, a ray, or an infinite line:
 * - If `lineType` is 1, it checks if the point lies on the **line segment** from A to B.
 * - If `lineType` is 2, it checks if the point lies on the **ray** starting at A and going through B.
 * - If `lineType` is 3, it checks if the point lies on the **infinite line** defined by A and B.
 *
 * @param[in] A The starting point of the line (1×3 vector).
 * @param[in] B The ending point of the line (1×3 vector).
 * @param[in] point The point to check (1×3 vector).
 * @param[in] lineType Specifies the type of line:
 *                     1 = Line segment,
 *                     2 = Ray extending from point A,
 *                     3 = Infinite line. (default: 1)
 * @return `true` if the point lies on the line based on `lineType`, `false` otherwise.
 */
bool isPointOnLine3D(
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    const Eigen::RowVector3d& point,
    int lineType
) {
    if (lineType < 1 || lineType > 3) {
        throw std::invalid_argument("Invalid lineType: must be 1 (segment), 2 (ray), or 3 (line)");
    }
    const Eigen::RowVector3d AB = B - A;
    const Eigen::RowVector3d AP = point - A;
    const double lengthSquaredAB = AB.squaredNorm();

    if (lengthSquaredAB < highEpsilon) return AP.squaredNorm() < highEpsilon;

    const double crossSquared = AB.cross(AP).squaredNorm();

    if (crossSquared > highEpsilon * lengthSquaredAB) return false;

    switch (lineType) {
    case 1: {
        const double AP_dot_AB = AP.dot(AB);
        return (AP_dot_AB >= -epsilon) && (AP_dot_AB <= lengthSquaredAB + epsilon);
    }
    case 2: {
        const double AP_dot_AB = AP.dot(AB);
        return AP_dot_AB >= -epsilon;
    }
    case 3:
        return true;
    }
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Determines whether a point lies inside a convex polygon in 2D.
 *
 * This function checks if a given point lies inside a 2D convex polygon or on its boundary.
 * If the `checkBoundary` parameter is set to true, the point is considered inside the polygon
 * even if it lies on the boundary.
 *
 * @param[in] polygonVertices The vertices of the convex polygon (N x 2 matrix, where N is the number of vertices).
 * @param[in] point The point to check (1 x 2 vector).
 * @param[in] checkBoundary If true, points on the boundary of the polygon are considered inside. (default: true)
 * @param[in] needClosePolygon  If true, ensures the polygon is closed by adding the first vertex to the end if not already closed. (default: true)
 *
 * @return `true` if the point lies inside the convex polygon or on the boundary (depending on `checkBoundary`), `false` otherwise.
 *
 */
bool isPointInPolygon2D(
    const Eigen::MatrixX2d& polygonVertices,
    const Eigen::RowVector2d& point,
    bool checkBoundary,
    bool needClosePolygon
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

    Eigen::Index numVertices = polygonPtr->rows() - 1;
    if (checkBoundary) {
        for (Eigen::Index i = 0; i < numVertices; ++i) {
            const auto& p1 = polygonPtr->row(i);
            const auto& p2 = polygonPtr->row(i + 1);
            if (isPointOnLine2D(p1, p2, point, 1)) {
                return true;
            }
        }
    }
    bool inside = false;
    for (Eigen::Index i = 0; i < numVertices; ++i) {
        const auto& p1 = polygonPtr->row(i);
        const auto& p2 = polygonPtr->row(i + 1);

        double y1 = p1(1), y2 = p2(1);
        double y_diff = y2 - y1;
        if (std::abs(y_diff) < epsilon) continue;

        double minY = std::min(y1, y2);
        double maxY = std::max(y1, y2);
        if (point(1) <= minY || point(1) > maxY) continue;

        double xIntersect = p1(0) + (point(1) - y1) * (p2(0) - p1(0)) / y_diff;
        if (point(0) <= xIntersect) inside = !inside;
    }
    return inside;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Determines whether a point lies inside a convex polygon in 3D.
 *
 * This function checks if a given point lies inside a 3D convex polygon or on its boundary.
 * If the `checkBoundary` parameter is set to true, the point is considered inside the polygon
 * even if it lies on the boundary.
 *
 * @param[in] polygonVertices The vertices of the convex polygon in 3D (N x 3 matrix, where N is the number of vertices).
 * @param[in] point The point to check (1 x 3 vector).
 * @param[in] checkBoundary If true, points on the boundary of the polygon are considered inside. (default: true)
 * @param[in] needClosePolygon If true, the polygon vertices will be sorted to form a closed polygon. (default: true)
 *
 * @return `true` if the point lies inside the convex polygon or on the boundary (depending on `checkBoundary`), `false` otherwise.
 *
 */
bool isPointInPolygon3D(
    const Eigen::MatrixX3d& polygonVertices,
    const Eigen::RowVector3d& point,
    bool checkBoundary,
    bool needClosePolygon
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
    Eigen::RowVector3d faceNormal = computePlaneNormal(polygonVertices).first;
    const Eigen::RowVector3d deltaPoint = point - refPoint;
    double pointPlaneDist = std::abs(faceNormal.dot(deltaPoint));
    if (pointPlaneDist > epsilon) return false;

    Eigen::Matrix<double, 3, 2> T = computeProjectionMatrix(faceNormal);

    Eigen::RowVector2d projectedPoint = deltaPoint * T;
    Eigen::MatrixX2d projectedVertices = (polygonPtr->rowwise() - refPoint) * T;

    return isPointInPolygon2D(projectedVertices, projectedPoint, checkBoundary, false);
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Determines whether a point lies inside a convex polyhedron.
 *
 * This function checks if a given point lies inside a convex polyhedron defined by its vertices.
 * If the `checkBoundary` parameter is set to true, the point is considered inside the polyhedron
 * even if it lies on the boundary.
 *
 * @param[in] polygonVertices The vertices of the convex polyhedron (N x 3 matrix, where N is the number of vertices).
 * @param[in] point The point to check (1 x 3 vector).
 * @param[in] checkBoundary If true, points on the boundary of the polyhedron are considered inside. (default: true)
 *
 * @return `true` if the point lies inside the polyhedron (or on the boundary if `checkBoundary` is true),
 *         `false` otherwise.
 *
 */
bool isPointInPolyhedron(
    const Eigen::MatrixX3d& polygonVertices,
    const Eigen::RowVector3d& point,
    bool checkBoundary
) {
    Mesh polyhedron_mesh = buildConvexHullMesh(polygonVertices);
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

//************************************************************************************************************************//
/**
 * @brief Checks if two line segments or rays intersect and calculates the intersection point(s).
 *
 * This function checks whether two 2D line segments or rays intersect. If they do, it calculates the intersection point(s).
 * If `intersection` is provided, all intersection points are calculated and stored in the `intersection` variable.
 * If `intersection` is not provided, the function checks if the lines intersect and returns `true` as soon as the first intersection is found.
 *
 * @param[in] A The starting point of line segment/ray AB (2D vector, x and y coordinates).
 * @param[in] B The ending point of line segment/ray AB (2D vector, x and y coordinates).
 * @param[in] C The starting point of line segment/ray CD (2D vector, x and y coordinates).
 * @param[in] D The ending point of line segment/ray CD (2D vector, x and y coordinates).
 * @param[out] intersection If provided, the intersection points are stored in this variable.
 * @param[in] lineType Type of the lines:
 *                     - 1: Line segments,
 *                     - 2: Rays,
 *                     - 3: Infinite lines.
 * @return `true` if the lines intersect at least once; `false` if the lines do not intersect.
 *         If `intersection` is provided, all intersection points are stored in `intersection`.
 *
 */
bool isLinesIntersection2D(
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    const Eigen::RowVector2d& C,
    const Eigen::RowVector2d& D,
    Eigen::RowVector2d& intersection,
    int lineType
) {
    if (lineType < 1 || lineType > 3) {
        throw std::invalid_argument("Invalid lineType: must be 1 (segment), 2 (ray), or 3 (line)");
    }

    const Eigen::RowVector2d AB = B - A;
    const Eigen::RowVector2d CD = D - C;
    const Eigen::RowVector2d AC = C - A;

    const double denom = AB.cross(CD);
    if (std::abs(denom) < epsilon) return false;

    const double t = AC.cross(CD) / denom;
    const double u = AC.cross(AB) / denom;

    switch (lineType) {
    case 1:
        if (t < -epsilon || t > epsilon_plus_1 || u < -epsilon || u > epsilon_plus_1) return false;
    case 2:
        if (t < -epsilon || u < -epsilon) return false;
    case 3:
        break;
    }

    intersection = A + t * AB;
    return true;
}

/**
 * @brief Checks if two line segments or rays intersect and returns `true` as soon as the first intersection is found.
 *
 * This function checks if two 2D line segments or rays intersect. It returns `true` as soon as the first intersection is found.
 * The intersection point is not calculated; the function simply checks for the presence of an intersection.
 *
 * @param[in] A The starting point of line segment/ray AB (2D vector, x and y coordinates).
 * @param[in] B The ending point of line segment/ray AB (2D vector, x and y coordinates).
 * @param[in] C The starting point of line segment/ray CD (2D vector, x and y coordinates).
 * @param[in] D The ending point of line segment/ray CD (2D vector, x and y coordinates).
 * @param[in] lineType Type of the lines:
 *                     - 1: Line segments,
 *                     - 2: Rays,
 *                     - 3: Infinite lines.
 * @return `true` if the lines intersect at least once; `false` if the lines do not intersect.
 *
 */
bool isLinesIntersection2D(
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    const Eigen::RowVector2d& C,
    const Eigen::RowVector2d& D,
    int lineType
) {
    Eigen::RowVector2d intersection;
    return isLinesIntersection2D(A, B, C, D, intersection, lineType);
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Computes the shortest distance between two lines (or line segments/rays) in 3D space.
 *
 * This function calculates the closest points between two 3D lines (or line segments/rays),
 * and returns the two closest points along with the distance between them. The lines or rays
 * are represented by two pairs of points A-B and C-D.
 *
 * @param[in] A The starting point of the first line/segment/ray (3D vector).
 * @param[in] B The ending point of the first line/segment/ray (3D vector).
 * @param[in] C The starting point of the second line/segment/ray (3D vector).
 * @param[in] D The ending point of the second line/segment/ray (3D vector).
 * @param[in] lineType Specifies the type of the lines:
 *                     - 1: Line segments
 *                     - 2: Rays
 *                     - 3: Infinite lines
 * @return A tuple containing:
 *         - The closest point on the first line/segment/ray.
 *         - The closest point on the second line/segment/ray.
 *         - The shortest distance between the two lines/segments/rays.
 *
 */
std::tuple<Eigen::RowVector3d, Eigen::RowVector3d, double> computeLinesDistance(
    const Eigen::RowVector3d& A, const Eigen::RowVector3d& B,
    const Eigen::RowVector3d& C, const Eigen::RowVector3d& D,
    int lineType
) {
    if (lineType < 1 || lineType > 3) {
        throw std::invalid_argument("Invalid lineType: must be 1 (segment), 2 (ray), or 3 (line)");
    }

    const auto clampValue = [lineType](double value) {
        switch (lineType) {
        case 1: return std::clamp(value, 0.0, 1.0);
        case 2: return std::max(value, 0.0);
        case 3: return value;
        }
        };

    const Eigen::RowVector3d AB = B - A;
    const Eigen::RowVector3d CD = D - C;
    const Eigen::RowVector3d AC = C - A;

    const double lengthAB = AB.squaredNorm();
    const double lengthCD = CD.squaredNorm();
    double t = 0.0, u = 0.0;
    const bool is_degenerate = (lengthAB < highEpsilon) || (lengthCD < highEpsilon);
    if (is_degenerate) {
        if (lengthAB > highEpsilon) {
            t = clampValue(AB.dot(AC) / lengthAB);
        }
        else if (lengthCD > highEpsilon) {
            u = clampValue(-CD.dot(AC) / lengthCD);
        }
    }
    else {
        const double AB_dot_CD = AB.dot(CD);
        const double denominator = lengthAB * lengthCD - AB_dot_CD * AB_dot_CD;
        if (std::abs(denominator) < highEpsilon) {// 平行或共线
            u = clampValue(-CD.dot(AC) / lengthCD);
            t = clampValue((u * AB_dot_CD + AB.dot(AC)) / lengthAB);
        }
        else {// 一般情况
            t = clampValue((AB.dot(AC) * lengthCD - CD.dot(AC) * AB_dot_CD) / denominator);
            u = clampValue((t * AB_dot_CD - CD.dot(AC)) / lengthCD);
        }
    }

    Eigen::RowVector3d closestOnAB = A + AB * t;
    Eigen::RowVector3d closestOnCD = C + CD * u;
    double distance = (closestOnAB - closestOnCD).norm();

    return { closestOnAB, closestOnCD, distance };
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if two 3D lines (or line segments/rays) intersect within a specified distance threshold.
 *
 * This function uses the `computeLinesDistance` function to calculate the shortest distance
 * between two lines or line segments/rays. If the distance is less than or equal to the
 * given threshold, the lines are considered intersecting.
 *
 * @param[in] A The starting point of the first line/segment/ray (3D vector).
 * @param[in] B The ending point of the first line/segment/ray (3D vector).
 * @param[in] C The starting point of the second line/segment/ray (3D vector).
 * @param[in] D The ending point of the second line/segment/ray (3D vector).
 * @param[out] output A tuple that stores the closest points on the two lines and the distance.
 * @param[in] threshold The maximum allowable distance for the lines to be considered intersecting.
 * @param[in] lineType Specifies the type of the lines:
 *                     - 1: Line segments
 *                     - 2: Rays
 *                     - 3: Infinite lines
 * @return `true` if the distance between the two lines is less than or equal to the threshold,
 *         `false` otherwise.
 */
bool isLinesIntersection3D(
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    const Eigen::RowVector3d& C,
    const Eigen::RowVector3d& D,
    std::tuple<Eigen::RowVector3d, Eigen::RowVector3d, double>& output,
    double threshold,
    int lineType
) {
    output = computeLinesDistance(
        A, B, C, D, lineType
    );
    return std::get<2>(output) <= threshold + epsilon;
}

/**
 * @brief Checks if two 3D lines (or line segments/rays) intersect within a specified distance threshold.
 *
 * This function calculates the shortest distance between two lines/segments/rays and checks if
 * the distance is less than or equal to the given threshold. It uses the `computeLinesDistance`
 * function for the calculation, but does not return the closest points.
 *
 * @param[in] A The starting point of the first line/segment/ray (3D vector).
 * @param[in] B The ending point of the first line/segment/ray (3D vector).
 * @param[in] C The starting point of the second line/segment/ray (3D vector).
 * @param[in] D The ending point of the second line/segment/ray (3D vector).
 * @param[in] threshold The maximum allowable distance for the lines to be considered intersecting.
 * @param[in] lineType Specifies the type of the lines:
 *                     - 1: Line segments
 *                     - 2: Rays
 *                     - 3: Infinite lines
 * @return `true` if the distance between the two lines is less than or equal to the threshold,
 *         `false` otherwise.
 */
bool isLinesIntersection3D(
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    const Eigen::RowVector3d& C,
    const Eigen::RowVector3d& D,
    double threshold,
    int lineType
) {
    return std::get<2>(computeLinesDistance(A, B, C, D, lineType)) <= threshold + epsilon;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if a 2D line intersects with a polygon and returns the intersection points.
 *
 * This function checks if a given line intersects with the edges of a polygon. It returns
 * the intersection points between the line and the polygon's edges and stores them in the
 * provided vector `intersections`. The polygon is represented by its vertices in a 2D matrix.
 * The line is defined by two points A and B. The function also allows sorting the vertices of
 * the polygon before checking for intersections.
 *
 * @param[in] polygonVertices A matrix representing the vertices of the polygon (2D points).
 * @param[in] A The starting point of the line (2D vector).
 * @param[in] B The ending point of the line (2D vector).
 * @param[out] intersections A vector that will store the intersection points between the line and the polygon.
 * @param[in] lineType Specifies the type of the line:
 *                     - 1: Line segment
 *                     - 2: Ray
 *                     - 3: Infinite line
 * @param[in] needClosePolygon A boolean flag to indicate if the polygon vertices should be sorted before checking for intersections.
 * @return `true` if there are intersection points, `false` otherwise.
 */
bool isLinePolygonIntersection2D(
    const Eigen::MatrixX2d& polygonVertices,
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    std::vector<Eigen::RowVector2d>& intersections,
    int lineType,
    bool needClosePolygon
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
        if (isLinesIntersection2D(A, B, p1, p2, intersection, lineType)) {
            intersections.emplace_back(intersection);
        }
    }
    return !intersections.empty();
}

/**
 * @brief Checks if a 2D line intersects with a polygon.
 *
 * This function checks if a given line intersects with any of the edges of a polygon. The polygon
 * is represented by its vertices in a 2D matrix, and the line is defined by two points A and B.
 * The function also allows sorting the vertices of the polygon before checking for intersections.
 *
 * @param[in] polygonVertices A matrix representing the vertices of the polygon (2D points).
 * @param[in] A The starting point of the line (2D vector).
 * @param[in] B The ending point of the line (2D vector).
 * @param[in] lineType Specifies the type of the line:
 *                     - 1: Line segment
 *                     - 2: Ray
 *                     - 3: Infinite line
 * @param[in] needClosePolygon A boolean flag to indicate if the polygon vertices should be sorted before checking for intersections.
 * @return `true` if the line intersects with the polygon, `false` otherwise.
 */
bool isLinePolygonIntersection2D(
    const Eigen::MatrixX2d& polygonVertices,
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    int lineType,
    bool needClosePolygon
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
        Eigen::RowVector2d intersection;
        if (isLinesIntersection2D(A, B, p1, p2, intersection, lineType)) {
            return true;
        }
    }
    return false;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if a 3D line intersects with a 3D polygon and returns the intersection points.
 *
 * This function checks if a line intersects with the edges of a polygon in 3D space. It first ensures
 * that the polygon and the line are coplanar. Then, depending on the specified line type and the
 * projection of the polygon, it calculates the intersection points and stores them in the provided
 * vector `intersections`. The polygon is represented by its vertices, and the line is defined by two points A and B.
 *
 * @param[in] polygonVertices A matrix representing the vertices of the polygon in 3D space.
 * @param[in] A The starting point of the line in 3D space (vector).
 * @param[in] B The ending point of the line in 3D space (vector).
 * @param[out] intersections A vector that will store the intersection points between the line and the polygon.
 * @param[in] lineType Specifies the type of the line:
 *                     - 1: Line segment
 *                     - 2: Ray
 *                     - 3: Infinite line
 * @param[in] checkBoundary A boolean flag to indicate whether to check boundary conditions for intersections.
 * @param[in] needClosePolygon A boolean flag to indicate if the polygon vertices should be sorted before checking for intersections.
 * @return `true` if there are intersection points, `false` otherwise.
 */
bool isLinePolygonIntersection3D(
    const Eigen::MatrixX3d& polygonVertices,
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    std::vector<Eigen::RowVector3d>& intersections,
    int lineType,
    bool checkBoundary,
    bool needClosePolygon
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
    const Eigen::RowVector3d lineDirection = B - A;
    double denominator = faceNormal.dot(lineDirection);
    const Eigen::RowVector3d deltaA = A - refPoint;

    if (checkBoundary && std::abs(denominator) < epsilon) {
        const Eigen::RowVector3d deltaB = B - refPoint;
        if (std::abs(faceNormal.dot(deltaA)) > epsilon || std::abs(faceNormal.dot(deltaB)) > epsilon) {
            return false; // 线段不与平面共面
        }
        const Eigen::MatrixX3d deltaVertices = polygonPtr->rowwise() - refPoint;
        Eigen::Matrix<double, 3, 2> T = computeProjectionMatrix(faceNormal);
        std::vector<Eigen::RowVector2d> intersections2D;
        if (isLinePolygonIntersection2D(deltaVertices * T, deltaA * T, deltaB * T, intersections2D, lineType, false)) {
            for (const auto& p2D : intersections2D) {
                Eigen::RowVector3d p3D = refPoint + p2D * T.transpose();
                intersections.emplace_back(p3D);
            }
            return true;
        }
        else {
            return false;
        }
    }

    const  double t = faceNormal.dot(-deltaA) / denominator;

    switch (lineType) {
    case 1:
        if (t < -epsilon || t > epsilon_plus_1) return false; // 修正范围
        break;
    case 2:
        if (t < -epsilon) return false;
        break;
    case 3:
        break;
    }

    Eigen::RowVector3d intersection = A + t * lineDirection;
    Eigen::Matrix<double, 3, 2> T = computeProjectionMatrix(faceNormal);
    Eigen::RowVector2d projectedPoint = intersection * T;
    Eigen::MatrixX2d projectedVertices = (*polygonPtr) * T;
    if (isPointInPolygon2D(projectedVertices, projectedPoint, checkBoundary, false)) {
        intersections.emplace_back(intersection);
        return true;
    }
    return false;
}

/**
 * @brief Checks if a 3D line intersects with a 3D polygon and returns a boolean indicating the intersection status.
 *
 * This function is a simpler version of `isLinePolygonIntersection3D` that only returns a boolean value.
 * It calls the more complex version of the function and passes an empty vector to store the intersection points.
 *
 * @param[in] polygonVertices A matrix representing the vertices of the polygon in 3D space.
 * @param[in] A The starting point of the line in 3D space (vector).
 * @param[in] B The ending point of the line in 3D space (vector).
 * @param[in] lineType Specifies the type of the line:
 *                     - 1: Line segment
 *                     - 2: Ray
 *                     - 3: Infinite line
 * @param[in] checkBoundary A boolean flag to indicate whether to check boundary conditions for intersections.
 * @param[in] needClosePolygon A boolean flag to indicate if the polygon vertices should be sorted before checking for intersections.
 * @return `true` if the line intersects with the polygon, `false` otherwise.
 */
bool isLinePolygonIntersection3D(
    const Eigen::MatrixX3d& polygonVertices,
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    int lineType,
    bool checkBoundary,
    bool needClosePolygon
) {
    std::vector<Eigen::RowVector3d> intersections;
    return isLinePolygonIntersection3D(polygonVertices, A, B, intersections, lineType, checkBoundary, needClosePolygon);
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Checks if a 3D line intersects with a 3D polyhedron and stores the intersection points.
 *
 * This function checks if a line (or ray or segment) intersects with the faces of a polyhedron represented by a set of vertices.
 * If intersections are found, they are stored in the `intersections` vector.
 *
 * @param[in] polyhedronVertices A matrix representing the vertices of the polyhedron in 3D space.
 * @param[in] A The starting point of the line in 3D space (vector).
 * @param[in] B The ending point of the line in 3D space (vector).
 * @param[out] intersections A vector to store the intersection points in 3D space (if any).
 * @param[in] lineType Specifies the type of the line:
 *                     - 1: Line segment
 *                     - 2: Ray
 *                     - 3: Infinite line
 * @return `true` if the line intersects with the polyhedron, `false` otherwise.
 */
bool isLinePolyhedronIntersection(
    const Eigen::MatrixX3d& polyhedronVertices,
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    std::vector<Eigen::RowVector3d>& intersections,
    int lineType
) {
    if (lineType < 1 || lineType > 3) {
        throw std::invalid_argument("Invalid lineType: must be 1 (segment), 2 (ray), or 3 (line)");
    }

    const Eigen::RowVector3d lineDirection = B - A;
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
        double denominator = faceNormal.dot(lineDirection);
        if (std::abs(denominator) < epsilon) continue;

        const double D = faceNormal.dot(eigenPoly.row(0));
        const double t = (D - faceNormal.dot(A)) / denominator;

        switch (lineType) {
        case 1:
            if (t < -epsilon || t > epsilon_plus_1) continue;
            break;
        case 2:
            if (t < -epsilon) continue;
            break;
        case 3:
            break;
        }

        Eigen::RowVector3d intersection = A + t * lineDirection;
        Eigen::Matrix<double, 3, 2> T = computeProjectionMatrix(faceNormal);

        Eigen::RowVector2d projectedPoint = intersection * T;
        Eigen::MatrixX2d projectedVertices = eigenPoly * T;

        if (isPointInPolygon2D(projectedVertices, projectedPoint, false, false)) {
            intersections.emplace_back(intersection);
        }
    }
    return !intersections.empty();
}

/**
 * @brief Checks if a 3D line intersects with a 3D polyhedron and returns a boolean indicating the intersection status.
 *
 * This function is a simpler version of `isLinePolyhedronIntersection` that only returns a boolean value.
 * It checks if the line (or ray or segment) intersects any face of the polyhedron represented by a set of vertices.
 *
 * @param[in] polyhedronVertices A matrix representing the vertices of the polyhedron in 3D space.
 * @param[in] A The starting point of the line in 3D space (vector).
 * @param[in] B The ending point of the line in 3D space (vector).
 * @param[in] lineType Specifies the type of the line:
 *                     - 1: Line segment
 *                     - 2: Ray
 *                     - 3: Infinite line
 * @return `true` if the line intersects with the polyhedron, `false` otherwise.
 */
bool isLinePolyhedronIntersection(
    const Eigen::MatrixX3d& polyhedronVertices,
    const Eigen::RowVector3d& A,
    const Eigen::RowVector3d& B,
    int lineType
) {
    if (lineType < 1 || lineType > 3) {
        throw std::invalid_argument("Invalid lineType: must be 1 (segment), 2 (ray), or 3 (line)");
    }

    const Eigen::RowVector3d lineDirection = B - A;
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
        double denominator = faceNormal.dot(lineDirection);
        if (std::abs(denominator) < epsilon) continue;

        const double D = faceNormal.dot(eigenPoly.row(0));
        const double t = (D - faceNormal.dot(A)) / denominator;

        switch (lineType) {
        case 1:
            if (t < -epsilon || t > epsilon_plus_1) continue;
            break;
        case 2:
            if (t < -epsilon) continue;
            break;
        case 3:
            break;
        }

        Eigen::RowVector3d intersection = A + t * lineDirection;
        Eigen::Matrix<double, 3, 2> T = computeProjectionMatrix(faceNormal);

        Eigen::RowVector2d projectedPoint = intersection * T;
        Eigen::MatrixX2d projectedVertices = eigenPoly * T;
        if (isPointInPolygon2D(projectedVertices, projectedPoint, false, false)) {
            return true;
        }
    }
    return false;
}

//************************************************************************************************************************//

//************************************************************************************************************************//
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
 * The function uses parallel processing for large datasets to improve efficiency.
 *
 * @param[in] polygonVertices A matrix representing the vertices of the polygon in 2D space.
 * @param[in] points A vector of 2D points to be filtered.
 * @param[in] removeBoundary A boolean flag to control whether boundary points are removed.
 * @param[in] removeInside A boolean flag to control whether inside points are removed.
 * @param[in] needClosePolygon A boolean flag to control whether the polygon vertices should be sorted before filtering.
 * @return A vector containing the filtered points.
 */
std::vector<Eigen::RowVector2d> filterPointsByPolygon(
    const Eigen::MatrixX2d& polygonVertices,
    const std::vector<Eigen::RowVector2d>& points,
    bool removeBoundary,
    bool removeInside,
    bool needClosePolygon
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

    std::function<bool(const Eigen::RowVector2d&)> filterCondition;
    if (removeInside) {
        filterCondition = [&](const Eigen::RowVector2d& p) {
            return !isPointInPolygon2D(*polygonPtr, p, removeBoundary, false);
            };
    }
    else {
        filterCondition = [&](const Eigen::RowVector2d& p) {
            return isPointInPolygon2D(*polygonPtr, p, removeBoundary, false);
            };
    }

    size_t numPoints = points.size();
    std::vector<Eigen::RowVector2d> filteredPoints;
    filteredPoints.reserve(numPoints);
    if (numPoints < 10000) {
        filteredPoints.reserve(numPoints);
        for (const auto& p : points) {
            if (filterCondition(p)) {
                filteredPoints.emplace_back(p);
            }
        }
    }
    else {
        size_t hardwareThreads = std::thread::hardware_concurrency();
        size_t minChunkSize = 100;
        size_t maxChunkSize = 10000;

        size_t dynamicChunkSize = std::clamp(
            numPoints / (4 * hardwareThreads),
            minChunkSize,
            maxChunkSize
        );

        size_t numChunks = (numPoints + dynamicChunkSize - 1) / dynamicChunkSize;
        std::vector<std::vector<Eigen::RowVector2d>> chunkResults(numChunks);

        std::vector<size_t> chunkIndices(numChunks);
        std::iota(chunkIndices.begin(), chunkIndices.end(), 0);

        std::for_each(std::execution::par, chunkIndices.begin(), chunkIndices.end(), [&](size_t chunkIdx) {
            size_t start = chunkIdx * dynamicChunkSize;
            size_t end = std::min(start + dynamicChunkSize, numPoints);
            auto& local = chunkResults[chunkIdx];
            local.reserve(end - start);
            for (size_t i = start; i < end; ++i) {
                if (filterCondition(points[i])) {
                    local.emplace_back(points[i]);
                }
            }
            });

        size_t total = 0;
        for (const auto& vec : chunkResults) total += vec.size();
        filteredPoints.reserve(total);

        for (auto& vec : chunkResults) {
            filteredPoints.insert(filteredPoints.end(),
                std::make_move_iterator(vec.begin()),
                std::make_move_iterator(vec.end()));
        }
    }
    filteredPoints.shrink_to_fit();
    return filteredPoints;
}
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
 * The function utilizes the CGAL library's `Point_inside` detection to determine if points are inside or outside the polyhedron.
 * Parallel processing is used for large datasets to improve performance.
 *
 * @param[in] polyhedronVertices A matrix representing the vertices of the polyhedron in 3D space.
 * @param[in] points A vector of 3D points to be filtered.
 * @param[in] removeBoundary A boolean flag to control whether boundary points are removed.
 * @param[in] removeInside A boolean flag to control whether inside points are removed.
 * @return A vector containing the filtered points.
 */
std::vector<Eigen::RowVector3d> filterPointsByPolyhedron(
    const Eigen::MatrixX3d& polyhedronVertices,
    const std::vector<Eigen::RowVector3d>& points,
    bool removeBoundary,
    bool removeInside
) {
    Mesh polyhedronMesh = buildConvexHullMesh(polyhedronVertices);
    auto meshFaceIterators = faces(polyhedronMesh);
    Tree tree(
        std::make_move_iterator(meshFaceIterators.first),
        std::make_move_iterator(meshFaceIterators.second),
        polyhedronMesh);

    Point_inside insideDetecter(std::move(tree));
    std::function<bool(const Eigen::RowVector3d&)> filterCondition;
    if (!removeInside && !removeBoundary) {
        // 保留内部和边界，过滤外部：只需要判断是否非外部
        filterCondition = [&](const Eigen::RowVector3d& p) {
            return insideDetecter(Point_3(p.x(), p.y(), p.z())) != CGAL::ON_UNBOUNDED_SIDE;
            };
    }
    else if (!removeInside && removeBoundary) {
        // 保留内部，过滤边界和外部：严格判断是否为内部
        filterCondition = [&](const Eigen::RowVector3d& p) {
            return insideDetecter(Point_3(p.x(), p.y(), p.z())) == CGAL::ON_BOUNDED_SIDE;
            };
    }
    else if (removeInside && !removeBoundary) {
        // 保留外部和边界，过滤内部：判断是否非内部
        filterCondition = [&](const Eigen::RowVector3d& p) {
            return insideDetecter(Point_3(p.x(), p.y(), p.z())) != CGAL::ON_BOUNDED_SIDE;
            };
    }
    else {
        // 保留外部，过滤内部和边界：严格判断是否为外部
        filterCondition = [&](const Eigen::RowVector3d& p) {
            return insideDetecter(Point_3(p.x(), p.y(), p.z())) == CGAL::ON_UNBOUNDED_SIDE;
            };
    }

    std::vector<Eigen::RowVector3d> filteredPoints;
    size_t numPoints = points.size();
    filteredPoints.reserve(numPoints);
    if (numPoints < 10000) {
        for (const auto& p : points) {
            if (filterCondition(p)) {
                filteredPoints.emplace_back(p);
            }
        }
    }
    else {
        size_t hardwareThreads = std::thread::hardware_concurrency();
        size_t minChunkSize = 100;
        size_t maxChunkSize = 10000;

        size_t dynamicChunkSize = std::clamp(
            numPoints / (4 * hardwareThreads),
            minChunkSize,
            maxChunkSize
        );

        size_t numChunks = (numPoints + dynamicChunkSize - 1) / dynamicChunkSize;
        std::vector<std::vector<Eigen::RowVector3d>> chunkResults(numChunks);

        std::vector<size_t> chunkIndices(numChunks);
        std::iota(chunkIndices.begin(), chunkIndices.end(), 0);

        std::for_each(std::execution::par, chunkIndices.begin(), chunkIndices.end(), [&](size_t chunkIdx) {
            size_t start = chunkIdx * dynamicChunkSize;
            size_t end = std::min(start + dynamicChunkSize, numPoints);
            auto& local = chunkResults[chunkIdx];
            local.reserve(end - start);
            for (size_t i = start; i < end; ++i) {
                if (filterCondition(points[i])) {
                    local.emplace_back(points[i]);
                }
            }
            });

        size_t total = 0;
        for (const auto& vec : chunkResults) total += vec.size();
        filteredPoints.reserve(total);

        for (auto& vec : chunkResults) {
            filteredPoints.insert(filteredPoints.end(),
                std::make_move_iterator(vec.begin()),
                std::make_move_iterator(vec.end()));
        }
    }
    filteredPoints.shrink_to_fit();
    return filteredPoints;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
static Eigen::MatrixX2d computeConvexHull2D(const Eigen::MatrixX2d& points) {
    auto cross = [](
        const Eigen::RowVector2d& a,
        const Eigen::RowVector2d& b,
        const Eigen::RowVector2d& c) -> double {
            return (b.x() - a.x()) * (c.y() - a.y()) -
                (b.y() - a.y()) * (c.x() - a.x());
        };

    size_t numPoints = points.rows();
    if (numPoints < 3) return points;

    std::vector<Eigen::RowVector2d> sortedPoints(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        sortedPoints[i] = points.row(i);
    }

    std::sort(sortedPoints.begin(), sortedPoints.end(),
        [](const Eigen::RowVector2d& a, const Eigen::RowVector2d& b) {
            return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y());
        });

    auto last = std::unique(sortedPoints.begin(), sortedPoints.end(),
        [](const Eigen::RowVector2d& a, const Eigen::RowVector2d& b) {
            return a.isApprox(b);
        });

    sortedPoints.erase(last, sortedPoints.end());

    std::vector<Eigen::RowVector2d> hull;
    hull.reserve(numPoints);
    for (size_t i = 0; i < sortedPoints.size(); ++i) {
        while (hull.size() >= 2 &&
            cross(hull[hull.size() - 2], hull.back(), sortedPoints[i]) <= 0) {
            hull.pop_back();
        }
        hull.emplace_back(sortedPoints[i]);
    }

    // 构建上凸包
    size_t lowerLen = hull.size();
    for (int i = static_cast<int>(sortedPoints.size()) - 2; i >= 0; --i) {
        while (hull.size() > lowerLen &&
            cross(hull[hull.size() - 2], hull.back(), sortedPoints[i]) <= 0) {
            hull.pop_back();
        }
        hull.emplace_back(sortedPoints[i]);
    }

    if (!hull.empty() && hull.front() == hull.back()) {
        hull.pop_back();
    }

    if (hull.size() == 2 && hull[0] == hull[1]) {
        hull.resize(1);
    }

    Eigen::MatrixX2d hullVertices(hull.size(), 2);
    for (size_t i = 0; i < hull.size(); ++i) {
        hullVertices.row(i) = hull[i];
    }

    return hullVertices;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
static Eigen::MatrixX3d computeConvexHull3D(const Eigen::MatrixX3d& points) {
    Eigen::Index numPoints = points.rows();
    if (numPoints < 4) {
        return points;
    }

    std::vector<Point_3> cgalPoints;
    cgalPoints.reserve(numPoints);
    for (Eigen::Index i = 0; i < numPoints; ++i) {
        cgalPoints.emplace_back(points(i, 0), points(i, 1), points(i, 2));
    }

    Polyhedron_3 hull;
    CGAL::convex_hull_3(cgalPoints.begin(), cgalPoints.end(), hull);
    Eigen::MatrixX3d hull_vertices(hull.size_of_vertices(), 3);
    size_t idx = 0;
    for (auto vit = hull.vertices_begin(); vit != hull.vertices_end(); ++vit, ++idx) {
        const auto& p = vit->point();
        hull_vertices(idx, 0) = CGAL::to_double(p.x());
        hull_vertices(idx, 1) = CGAL::to_double(p.y());
        hull_vertices(idx, 2) = CGAL::to_double(p.z());
    }
    return hull_vertices;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Computes the Principal Component Analysis (PCA) rectangle for a set of 2D points.
 *
 * This function computes the PCA rectangle that best fits the given set of 2D points.
 * The process includes optionally using the convex hull of the points, calculating the covariance matrix,
 * computing the eigenvectors (principal axes), and determining the bounding rectangle in the PCA space.
 * The result is a rectangle defined by four corner points, along with the centroid (mean) of the rectangle.
 *
 * @param[in] points A matrix of 2D points (each row is a point).
 * @param[in] useConvexHull A boolean flag indicating whether to compute the convex hull of the points.
 * @return A pair consisting of:
 *         - A 4x2 matrix representing the corners of the PCA rectangle.
 *         - A 2D vector representing the centroid (mean) of the rectangle.
 */
std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computePCARect(
    const Eigen::MatrixX2d& points,
    bool useConvexHull
) {
    Eigen::MatrixX2d processedPoints = useConvexHull ?
        computeConvexHull2D(points) : points;

    Eigen::Matrix2d covMatrix = computeCovarianceMatrix(processedPoints);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigenSolver(covMatrix);
    const Eigen::Matrix2d& eigenvectors = eigenSolver.eigenvectors();
    const Eigen::Vector2d& axis1 = eigenvectors.col(0);
    const Eigen::Vector2d& axis2 = eigenvectors.col(1);

    Eigen::MatrixX2d proj = processedPoints * eigenvectors;
    Eigen::Vector2d minProj = proj.colwise().minCoeff();
    Eigen::Vector2d maxProj = proj.colwise().maxCoeff();

    Eigen::Matrix<double, 4, 2> corners;
    corners.row(0) = axis1 * minProj(0) + axis2 * minProj(1);
    corners.row(1) = axis1 * maxProj(0) + axis2 * minProj(1);
    corners.row(2) = axis1 * maxProj(0) + axis2 * maxProj(1);
    corners.row(3) = axis1 * minProj(0) + axis2 * maxProj(1);

    return { corners, corners.colwise().mean() };
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Computes the Principal Component Analysis (PCA) box for a set of 3D points.
 *
 * This function computes the PCA box (axis-aligned bounding box in PCA space) that best fits the given set of 3D points.
 * The process includes optionally using the convex hull of the points, calculating the covariance matrix,
 * computing the eigenvectors (principal axes), and determining the bounding box in the PCA space.
 * The result is a box defined by eight corner points, along with the centroid (mean) of the box.
 *
 * @param[in] points A matrix of 3D points (each row is a point).
 * @param[in] useConvexHull A boolean flag indicating whether to compute the convex hull of the points.
 * @return A pair consisting of:
 *         - A 8x3 matrix representing the corners of the PCA box.
 *         - A 3D vector representing the centroid (mean) of the box.
 */
std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computePCABox(
    const Eigen::MatrixX3d& points,
    bool useConvexHull
) {
    Eigen::MatrixX3d processedPoints = useConvexHull ?
        computeConvexHull3D(points) : points;

    Eigen::Matrix3d covMatrix = computeCovarianceMatrix(processedPoints);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(covMatrix);
    const Eigen::Matrix3d& eigenvectors = eigenSolver.eigenvectors();
    const Eigen::Vector3d& axis1 = eigenvectors.col(0);
    const Eigen::Vector3d& axis2 = eigenvectors.col(1);
    const Eigen::Vector3d& axis3 = eigenvectors.col(2);
    Eigen::MatrixX3d proj = processedPoints * eigenvectors;

    Eigen::Vector3d minProj = proj.colwise().minCoeff();
    Eigen::Vector3d maxProj = proj.colwise().maxCoeff();

    Eigen::Matrix<double, 8, 3> corners;
    corners.row(0) = axis1 * minProj.x() + axis2 * minProj.y() + axis3 * minProj.z();
    corners.row(1) = axis1 * maxProj.x() + axis2 * minProj.y() + axis3 * minProj.z();
    corners.row(2) = axis1 * maxProj.x() + axis2 * maxProj.y() + axis3 * minProj.z();
    corners.row(3) = axis1 * minProj.x() + axis2 * maxProj.y() + axis3 * minProj.z();
    corners.row(4) = axis1 * minProj.x() + axis2 * minProj.y() + axis3 * maxProj.z();
    corners.row(5) = axis1 * maxProj.x() + axis2 * minProj.y() + axis3 * maxProj.z();
    corners.row(6) = axis1 * maxProj.x() + axis2 * maxProj.y() + axis3 * maxProj.z();
    corners.row(7) = axis1 * minProj.x() + axis2 * maxProj.y() + axis3 * maxProj.z();

    return { corners, corners.colwise().mean() };
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Computes the minimum bounding rectangle (MBR) for a set of 2D points.
 *
 * This function computes the minimum bounding rectangle (MBR) for the given set of 2D points.
 * The process includes calculating the convex hull of the points, and then determining the MBR that tightly fits the convex hull.
 * The result is a rectangle defined by four corner points, along with the centroid (mean) of the rectangle.
 *
 * @param[in] points A matrix of 2D points (each row is a point).
 * @return A pair consisting of:
 *         - A 4x2 matrix representing the corners of the minimum bounding rectangle (MBR).
 *         - A 2D vector representing the centroid (mean) of the rectangle.
 */
std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computeMinBoundRect(
    const Eigen::MatrixX2d& points
) {
    Eigen::Index numPoints = points.rows();
    std::vector<Point_2> cgalPoints;
    cgalPoints.reserve(numPoints);
    for (Eigen::Index i = 0; i < numPoints; ++i) {
        cgalPoints.emplace_back(points(i, 0), points(i, 1));
    }

    std::vector<Point_2> convexHull;
    convexHull.reserve(numPoints);
    CGAL::convex_hull_2(cgalPoints.begin(), cgalPoints.end(), std::back_inserter(convexHull));

    std::vector<Point_2> rectPoints;
    rectPoints.reserve(4);
    CGAL::min_rectangle_2(
        convexHull.begin(),
        convexHull.end(),
        std::back_inserter(rectPoints)
    );

    Eigen::Matrix<double, 4, 2> corners;
    for (size_t i = 0; i < 4; ++i) {
        corners.row(i) << CGAL::to_double(rectPoints[i].x()),
            CGAL::to_double(rectPoints[i].y());
    }

    return { corners, corners.colwise().mean() };
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Computes the minimum bounding box (MBB) for a set of 3D points.
 *
 * This function computes the minimum bounding box (MBB) for the given set of 3D points.
 * The process includes calculating the convex hull of the points, and then determining the oriented bounding box (OBB)
 * that tightly fits the convex hull in 3D space. The result is a box defined by eight corner points,
 * along with the centroid (mean) of the box.
 *
 * @param[in] points A matrix of 3D points (each row is a point).
 * @return A pair consisting of:
 *         - A 8x3 matrix representing the corners of the minimum bounding box (MBB).
 *         - A 3D vector representing the centroid (mean) of the box.
 */
std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computeMinBoundBox(
    const Eigen::MatrixX3d& points
) {
    Eigen::Index numPoints = points.rows();
    std::vector<Point_3> cgalPoints;
    cgalPoints.reserve(numPoints);
    for (Eigen::Index i = 0; i < numPoints; ++i) {
        cgalPoints.emplace_back(points(i, 0), points(i, 1), points(i, 2));
    }

    std::array<Point_3, 8> obbVertices;
    CGAL::oriented_bounding_box(cgalPoints, obbVertices,
        CGAL::parameters::use_convex_hull(true)
    );

    Eigen::Matrix<double, 8, 3> corners;
    for (size_t i = 0; i < 8; ++i) {
        corners.row(i) << CGAL::to_double(obbVertices[reorderIndexForCGAL[i]].x()),
            CGAL::to_double(obbVertices[reorderIndexForCGAL[i]].y()),
            CGAL::to_double(obbVertices[reorderIndexForCGAL[i]].z());
    }
    return { corners, corners.colwise().mean() };
}
//************************************************************************************************************************//

//************************************************************************************************************************//
/**
 * @brief Converts 3D world coordinates to 2D camera image coordinates.
 *
 * This function transforms a set of 3D target vertices from world coordinates to 2D camera image coordinates.
 * The transformation uses the camera position, the direction of the camera (from the camera position to the camera direction point),
 * and the focal length (FL) of the camera. The process involves computing the rotation matrix to align the world coordinates
 * with the camera coordinate system and then projecting the points onto the camera's image plane.
 *
 * @param[in] targetVertices A matrix of 3D target vertices (each row represents a vertex).
 * @param[in] cameraDirPoint A 3D point representing the direction of the camera.
 * @param[in] cameraPosition A 3D point representing the position of the camera.
 * @param[in] FL The focal length of the camera.
 * @return A matrix of 2D image coordinates corresponding to the input 3D target vertices.
 */
Eigen::MatrixX2d WorldToCameraImageCoords(
    const Eigen::MatrixX3d& targetVertices,
    const Eigen::RowVector3d& cameraDirPoint,
    const Eigen::RowVector3d& cameraPosition,
    double FL
) {
    Eigen::RowVector3d zc = (cameraDirPoint - cameraPosition).normalized();
    Eigen::RowVector3d xc = zc.cross(Eigen::RowVector3d::UnitZ()).normalized();
    Eigen::RowVector3d yc = xc.cross(zc).normalized();

    // calculate R.T
    Eigen::Matrix<double, 3, 3> R_T;
    R_T.row(0) << xc(0), yc(0), zc(0);
    R_T.row(1) << xc(1), yc(1), zc(1);
    R_T.row(2) << xc(2), yc(2), zc(2);

    Eigen::MatrixX3d cameraCoords = (targetVertices.rowwise() - cameraPosition) * R_T;

    Eigen::ArrayXd Z_inv = FL / cameraCoords.col(2).array();
    Eigen::MatrixX2d imageCoords(targetVertices.rows(), 2);

    imageCoords.col(0) = (cameraCoords.col(0).array() * Z_inv);  // u = FL*X/Z
    imageCoords.col(1) = (cameraCoords.col(1).array() * Z_inv); // v = FL*Y/Z

    return imageCoords;
}


/**
 * @brief Converts 3D world coordinates to 2D camera image coordinates (using the centroid as the camera direction point).
 *
 * This function is a simplified version of the previous function where the camera direction point is assumed to be
 * the centroid (mean) of the target vertices. It uses the camera position and focal length (FL) to convert the target
 * vertices from world coordinates to 2D camera image coordinates.
 *
 * @param[in] targetVertices A matrix of 3D target vertices (each row represents a vertex).
 * @param[in] cameraPosition A 3D point representing the position of the camera.
 * @param[in] FL The focal length of the camera.
 * @return A matrix of 2D image coordinates corresponding to the input 3D target vertices.
 */
Eigen::MatrixX2d WorldToCameraImageCoords(
    const Eigen::MatrixX3d& targetVertices,
    const Eigen::RowVector3d& cameraPosition,
    double FL
) {
    return WorldToCameraImageCoords(targetVertices, targetVertices.colwise().mean(), cameraPosition, FL);
}
//************************************************************************************************************************//