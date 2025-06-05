#include "RasterizationTool.h"
#include "PointTool.h"
static const std::array<Eigen::RowVector2i, 4> offsets2d = { {
    Eigen::RowVector2i(0, 0),  
    Eigen::RowVector2i(1, 0),
    Eigen::RowVector2i(0, 1),
    Eigen::RowVector2i(1, 1)
} };

static const std::array<Eigen::RowVector3i, 8> offsets3d = { {
    Eigen::RowVector3i(0, 0, 0),
    Eigen::RowVector3i(1, 0, 0),
    Eigen::RowVector3i(0, 1, 0),
    Eigen::RowVector3i(1, 1, 0),
    Eigen::RowVector3i(0, 0, 1),
    Eigen::RowVector3i(1, 0, 1),
    Eigen::RowVector3i(0, 1, 1),
    Eigen::RowVector3i(1, 1, 1)
} };

static std::vector<Eigen::RowVector2i> _rasterizeBresenhamLine2D(
    const Eigen::RowVector2d& p1,
    const Eigen::RowVector2d& p2
) {
    Eigen::RowVector2i p1_i = p1.array().floor().cast<int>();
    Eigen::RowVector2i p2_i = p2.array().floor().cast<int>();

    int x0 = p1_i.x();
    int y0 = p1_i.y();
    int x1 = p2_i.x();
    int y1 = p2_i.y();

    int dx = abs(x1 - x0);
    int dy = abs(y1 - y0);

    int maxSetp = std::max({ dx, dy }) + 1;
    std::vector<Eigen::RowVector2i> lineCells;
    lineCells.reserve(maxSetp);
    lineCells.emplace_back(x0, y0);

    int sx = x0 < x1 ? 1 : -1;
    int sy = y0 < y1 ? 1 : -1;
    int err = dx - dy;

    while (x0 != x1 || y0 != y1) {
        int e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x0 += sx;
        }
        if (e2 < dx) {
            err += dx;
            y0 += sy;
        }
        lineCells.emplace_back(x0, y0);
    }

    return lineCells;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
RowVector2iSet rasterizeBresenhamLine2D(
    const Eigen::RowVector2d& p1,
    const Eigen::RowVector2d& p2
) {
    std::vector<Eigen::RowVector2i> outputCells = _rasterizeBresenhamLine2D(p1, p2);
    RowVector2iSet gridPoints;
    gridPoints.reserve(outputCells.size() * 2);
    for (const auto& point : outputCells) {
        gridPoints.emplace(point);
        gridPoints.emplace(point + offsets2d[1]);
        gridPoints.emplace(point + offsets2d[2]);
        gridPoints.emplace(point + offsets2d[3]);
    }
    return gridPoints;
}


RowVector2iSet rasterizeBresenhamLine2D(
    const Eigen::RowVector2d& p1,
    const Eigen::RowVector2d& p2,
    std::vector<Eigen::RowVector2i>& rawLinePoints
) {
    rawLinePoints = _rasterizeBresenhamLine2D(p1, p2);
    RowVector2iSet gridPoints;
    gridPoints.reserve(rawLinePoints.size() * 2);

    for (const auto& point : rawLinePoints) {
        gridPoints.emplace(point);
        gridPoints.emplace(point + offsets2d[1]);
        gridPoints.emplace(point + offsets2d[2]);
        gridPoints.emplace(point + offsets2d[3]);
    }
    return gridPoints;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
RowVector2iSet rasterizePolygon2D(
    const Eigen::MatrixX2d& polygonVertices,
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

    RowVector2iSet boundaryGridPoints;
    boundaryGridPoints.reserve(1000);
    for (Eigen::Index i = 0; i < numVertices; i++) {
        const Eigen::RowVector2d& p1 = polygonPtr->row(i);
        const Eigen::RowVector2d& p2 = polygonPtr->row(i + 1);

        std::vector<Eigen::RowVector2i> outputCells;
        rasterizeBresenhamLine2D(p1, p2, outputCells);
        for (const auto& point : outputCells) {
            boundaryGridPoints.emplace(point);
        }
    }

    Eigen::RowVector2i min_vec = (polygonPtr->colwise().minCoeff().array().floor()).cast<int>();
    Eigen::RowVector2i max_vec = (polygonPtr->colwise().maxCoeff().array().ceil()).cast<int>();
    int x_min = min_vec(0), x_max = max_vec(0);
    int y_min = min_vec(1), y_max = max_vec(1);
    for (int x = x_min; x < x_max; x++) {
        for (int y = y_min; y < y_max; y++) {
            Eigen::RowVector2d center(x + 0.5, y + 0.5);
            if (isPointInPolygon2D(*polygonPtr, center, true, false)) {
                boundaryGridPoints.emplace(x, y);
            }
        }
    }

    RowVector2iSet polygonGridPoints;
    polygonGridPoints.reserve(boundaryGridPoints.size()*2);
    for (const auto& point : boundaryGridPoints) {
        polygonGridPoints.emplace(point);
        polygonGridPoints.emplace(point + offsets2d[1]);
        polygonGridPoints.emplace(point + offsets2d[2]);
        polygonGridPoints.emplace(point + offsets2d[3]);
    }
    return polygonGridPoints;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
static std::vector<Eigen::RowVector3i> _rasterizeBresenhamLine3D(
    const Eigen::RowVector3d& p1,
    const Eigen::RowVector3d& p2
) {
    Eigen::RowVector3i p1_i = p1.array().floor().cast<int>();
    Eigen::RowVector3i p2_i = p2.array().floor().cast<int>();

    int x0 = p1_i.x();
    int y0 = p1_i.y();
    int z0 = p1_i.z();
    int x1 = p2_i.x();
    int y1 = p2_i.y();
    int z1 = p2_i.z();

    int dx = abs(x1 - x0);
    int dy = abs(y1 - y0);
    int dz = abs(z1 - z0);

    int maxSetp = std::max({ dx, dy, dz }) + 1;
    std::vector<Eigen::RowVector3i> rawLinePoints;
    rawLinePoints.reserve(maxSetp);
    rawLinePoints.emplace_back(x0, y0, z0);

    int xs = (x1 > x0) ? 1 : -1;
    int ys = (y1 > y0) ? 1 : -1;
    int zs = (z1 > z0) ? 1 : -1;

    if (dx >= dy && dx >= dz) {
        // X-dominant
        int err_1 = 2 * dy - dx;
        int err_2 = 2 * dz - dx;
        for (int i = 0; i < dx; ++i) {
            x0 += xs;
            if (err_1 >= 0) {
                y0 += ys;
                err_1 -= 2 * dx;
            }
            if (err_2 >= 0) {
                z0 += zs;
                err_2 -= 2 * dx;
            }
            err_1 += 2 * dy;
            err_2 += 2 * dz;
            rawLinePoints.emplace_back(x0, y0, z0);
        }
    }
    else if (dy >= dx && dy >= dz) {
        // Y-dominant
        int err_1 = 2 * dx - dy;
        int err_2 = 2 * dz - dy;
        for (int i = 0; i < dy; ++i) {
            y0 += ys;
            if (err_1 >= 0) {
                x0 += xs;
                err_1 -= 2 * dy;
            }
            if (err_2 >= 0) {
                z0 += zs;
                err_2 -= 2 * dy;
            }
            err_1 += 2 * dx;
            err_2 += 2 * dz;
            rawLinePoints.emplace_back(x0, y0, z0);
        }
    }
    else {
        // Z-dominant
        int err_1 = 2 * dx - dz;
        int err_2 = 2 * dy - dz;
        for (int i = 0; i < dz; ++i) {
            z0 += zs;
            if (err_1 >= 0) {
                x0 += xs;
                err_1 -= 2 * dz;
            }
            if (err_2 >= 0) {
                y0 += ys;
                err_2 -= 2 * dz;
            }
            err_1 += 2 * dx;
            err_2 += 2 * dy;
            rawLinePoints.emplace_back(x0, y0, z0);
        }
    }
    return rawLinePoints;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
RowVector3iSet rasterizeBresenhamLine3D(
    const Eigen::RowVector3d& p1,
    const Eigen::RowVector3d& p2
) {
    std::vector<Eigen::RowVector3i> rawLinePoints = _rasterizeBresenhamLine3D(p1, p2);
    RowVector3iSet gridPoints;
    gridPoints.reserve(rawLinePoints.size() * 4);
    for (const auto& point : rawLinePoints) {
        gridPoints.emplace(point);
        gridPoints.emplace(point + offsets3d[1]);
        gridPoints.emplace(point + offsets3d[2]);
        gridPoints.emplace(point + offsets3d[3]);
        gridPoints.emplace(point + offsets3d[4]);
        gridPoints.emplace(point + offsets3d[5]);
        gridPoints.emplace(point + offsets3d[6]);
        gridPoints.emplace(point + offsets3d[7]);
    }
    return gridPoints;
}


RowVector3iSet rasterizeBresenhamLine3D(
    const Eigen::RowVector3d& p1,
    const Eigen::RowVector3d& p2,
    std::vector<Eigen::RowVector3i>& rawLinePoints
) {
    rawLinePoints = _rasterizeBresenhamLine3D(p1, p2);
    RowVector3iSet gridPoints;
    gridPoints.reserve(rawLinePoints.size() * 4);
    for (const auto& point : rawLinePoints) {
        gridPoints.emplace(point);
        gridPoints.emplace(point + offsets3d[1]);
        gridPoints.emplace(point + offsets3d[2]);
        gridPoints.emplace(point + offsets3d[3]);
        gridPoints.emplace(point + offsets3d[4]);
        gridPoints.emplace(point + offsets3d[5]);
        gridPoints.emplace(point + offsets3d[6]);
        gridPoints.emplace(point + offsets3d[7]);
    }
    return gridPoints;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
RowVector3iSet rasterizePolygon3D(
    const Eigen::MatrixX3d& polygonVertices,
    bool needSortVertices
) {
    Eigen::MatrixX3d polygon;
    const Eigen::MatrixX3d* polygonPtr = nullptr;
    if (needSortVertices) {
        polygon = generateClosedPolygon(polygonVertices);
        polygonPtr = &polygon;
    }
    else {
        polygonPtr = &polygonVertices;
    }
    Eigen::Index numVertices = polygonPtr->rows() - 1;

    RowVector3iSet boundaryGridPoints;
    boundaryGridPoints.reserve(1000);
    for (Eigen::Index i = 0; i < numVertices; ++i) {
        const auto& p1 = polygonPtr->row(i);
        const auto& p2 = polygonPtr->row(i + 1);

        std::vector<Eigen::RowVector3i> rawLinePoints;
        rasterizeBresenhamLine3D(p1, p2, rawLinePoints);

        for (const auto& point : rawLinePoints) {
            boundaryGridPoints.emplace(point);
        }
    }

    const Eigen::RowVector3d& origin = polygonPtr->row(0);
    const auto& [normal, _] = computePlaneNormal(polygonVertices);
    Eigen::Matrix<double, 3, 2> T = compute3Dto2DTransformMatrix(normal);
    Eigen::MatrixX2d polygon2D = (polygonPtr->rowwise() - origin) * T;

    Eigen::RowVector2i min_vec = (polygon2D.colwise().minCoeff().array().floor()).cast<int>();
    Eigen::RowVector2i max_vec = (polygon2D.colwise().maxCoeff().array().ceil()).cast<int>();
    int u_min = min_vec(0), u_max = max_vec(0);
    int v_min = min_vec(1), v_max = max_vec(1);

    for (int u = u_min; u < u_max; ++u) {
        for (int v = v_min; v < v_max; ++v) {
            Eigen::RowVector2d center2D(u + 0.5, v + 0.5);
            if (isPointInPolygon2D(polygon2D, center2D, true, false)) {
                Eigen::RowVector3i grid3D = (origin + center2D * T.transpose()).array().floor().cast<int>();
                boundaryGridPoints.emplace(grid3D);
            }
        }
    }

    RowVector3iSet polygonGridPoints;
    polygonGridPoints.reserve(boundaryGridPoints.size()*2);
    for (const auto& point : boundaryGridPoints) {
        polygonGridPoints.emplace(point);
        polygonGridPoints.emplace(point + offsets3d[1]);
        polygonGridPoints.emplace(point + offsets3d[2]);
        polygonGridPoints.emplace(point + offsets3d[3]);
        polygonGridPoints.emplace(point + offsets3d[4]);
        polygonGridPoints.emplace(point + offsets3d[5]);
        polygonGridPoints.emplace(point + offsets3d[6]);
        polygonGridPoints.emplace(point + offsets3d[7]);
    }

    return polygonGridPoints;
}
//************************************************************************************************************************//

//************************************************************************************************************************//

RowVector3iSet rasterizePolyhedron(
    const Eigen::MatrixX3d& polyhedronVertices
) {
    RowVector3iSet polyhedronGridPoints;
    Mesh poly = buildConvexHullMesh(polyhedronVertices);
    for (auto face : poly.faces()){
        auto he = poly.halfedge(face);
        auto start = he;

        Eigen::Matrix<double, 4, 3> eigenMatrix;
        Eigen::Index i = 0;
        do {
            const Point_3 point = poly.point(poly.target(he));
            eigenMatrix.row(i) << point.x(), point.y(), point.z();
            he = poly.next(he);
            i++;
        } while (he != start);
        eigenMatrix.row(3) = eigenMatrix.row(0);

        RowVector3iSet polygonGridPoints = rasterizePolygon3D(eigenMatrix, false);
        polyhedronGridPoints.merge(std::move(polygonGridPoints));
    }
    return polyhedronGridPoints;
}


RowVector3iSet rasterizePolyhedron(
    const Mesh& poly
) {
    RowVector3iSet polyhedronGridPoints;
    for (auto face : poly.faces()) {
        auto he = poly.halfedge(face);
        auto start = he;

        Eigen::Matrix<double, 4, 3> eigenMatrix;
        Eigen::Index i = 0;
        do {
            const Point_3 point = poly.point(poly.target(he));
            eigenMatrix.row(i) << point.x(), point.y(), point.z();
            he = poly.next(he);
            i++;
        } while (he != start);
        eigenMatrix.row(3) = eigenMatrix.row(0);
        //std::cout << "--------------------------" << std::endl;
        //std::cout << eigenMatrix.format(Eigen::FullPrecision) << std::endl;
        //std::cout << "--------------------------" << std::endl;
        RowVector3iSet polygonGridPoints = rasterizePolygon3D(eigenMatrix, false);
        polyhedronGridPoints.merge(std::move(polygonGridPoints));
    }
    return polyhedronGridPoints;
}
