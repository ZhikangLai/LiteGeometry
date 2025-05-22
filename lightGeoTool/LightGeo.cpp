#include "lightGeo.h"
#include <execution>
#include <numeric> 

static constexpr std::array<int, 8> reorderIndexForCGAL = { 0,1,2,3,5,6,7,4 };

static Eigen::MatrixXd computeCovarianceMatrix(
    const Eigen::MatrixXd& data,
    bool unbiased = true
) {
    Eigen::MatrixXd centered = data.rowwise() - data.colwise().mean();
    double divisor = unbiased ? static_cast<double>(data.rows() - 1) : static_cast<double>(data.rows());
    return (centered.transpose() * centered) / divisor;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
std::vector<Eigen::RowVector2d> filterPointsByPolygon(
    const Eigen::MatrixX2d& polygonVertices,
    const std::vector<Eigen::RowVector2d>& points,
    bool removeBoundary,
    bool removeInside,
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

    Mesh hull;
    CGAL::convex_hull_3(cgalPoints.begin(), cgalPoints.end(), hull);
    Eigen::MatrixX3d hull_vertices(hull.num_vertices(), 3);
    Eigen::Index idx = 0;

    for (auto v : hull.vertices()) {
        const auto& p = hull.point(v);
        hull_vertices(idx, 0) = CGAL::to_double(p.x());
        hull_vertices(idx, 1) = CGAL::to_double(p.y());
        hull_vertices(idx, 2) = CGAL::to_double(p.z());
        idx++;
    }

    return hull_vertices;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computePCAOBB2D(
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

std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computePCAOBB2D(
    const std::vector<Eigen::RowVector2d>& points,
    bool useConvexHull
) {
    Eigen::Index numPoints = points.size();
    Eigen::MatrixX2d _points(numPoints, 2);
    for (Eigen::Index i = 0; i < numPoints; ++i) {
        _points.row(i) = points[i];
    }
    return computePCAOBB2D(_points, useConvexHull);
}
//************************************************************************************************************************//

//************************************************************************************************************************//
std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computePCAOBB3D(
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

std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computePCAOBB3D(
    const std::vector<Eigen::RowVector3d>& points,
    bool useConvexHull
) {
    Eigen::Index numPoints = points.size();
    Eigen::MatrixX3d _points(numPoints, 3);
    for (Eigen::Index i = 0; i < numPoints; ++i) {
        _points.row(i) = points[i];
    }
    return computePCAOBB3D(_points, useConvexHull);
}
//************************************************************************************************************************//

//************************************************************************************************************************//
std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computeMinOBB2D(
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

std::pair<Eigen::Matrix<double, 4, 2>, Eigen::RowVector2d> computeMinOBB2D(
    const std::vector<Eigen::RowVector2d>& points
) {
    size_t numPoints = points.size();
    std::vector<Point_2> cgalPoints;
    cgalPoints.reserve(numPoints);
    for (const auto& point : points) {
        cgalPoints.emplace_back(point.x(), point.y());
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

std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computeMinOBB3D(
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


std::pair<Eigen::Matrix<double, 8, 3>, Eigen::RowVector3d> computeMinOBB3D(
    const std::vector<Eigen::RowVector3d>& points
) {
    size_t numPoints = points.size();
    std::vector<Point_3> cgalPoints;
    cgalPoints.reserve(numPoints);
    for (const auto& point : points) {
        cgalPoints.emplace_back(point.x(), point.y(), point.z());
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

Eigen::MatrixX2d WorldToCameraImageCoords(
    const Eigen::MatrixX3d& targetVertices,
    const Eigen::RowVector3d& cameraPosition,
    double FL
) {
    return WorldToCameraImageCoords(targetVertices, targetVertices.colwise().mean(), cameraPosition, FL);
}
//************************************************************************************************************************//