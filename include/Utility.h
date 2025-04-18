#ifndef UTILITY_H
#define UTILITY_H

#define EIGEN_USE_BLAS
#include <iostream>
#include <Eigen/Dense>
#include <unordered_set>
#include <unordered_map>
#include <algorithm> 
#include <string>
#include <vector>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Min_quadrilateral_2.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Point_2 Point_2;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> Tree;
typedef CGAL::Side_of_triangle_mesh<Mesh, K> Point_inside;
namespace PMP = CGAL::Polygon_mesh_processing;

namespace fs = std::filesystem;
static constexpr double inf = std::numeric_limits<double>::infinity();
static constexpr double epsilon = 1e-6;
static constexpr double highEpsilon = 1e-10;
static constexpr double epsilon_plus_1 = 1 + 1e-6;
static constexpr double PI = 3.141592653589793;
static constexpr double TWO_PI = 6.283185307179586;
static constexpr double RAD_TO_DEG = 57.2957786666617; // 180/pi

struct Segment2D {
    Eigen::RowVector2d P1, P2;
    Segment2D(const Eigen::RowVector2d& p1, const Eigen::RowVector2d& p2) : P1(p1), P2(p2) {}
};

struct Ray2D {
    Eigen::RowVector2d P1, P2;

    Ray2D(const Eigen::RowVector2d& p1, const Eigen::RowVector2d& p2) : P1(p1), P2(p2) {}
};

struct Line2D {
    Eigen::RowVector2d P1, P2;
    Line2D(const Eigen::RowVector2d& p1, const Eigen::RowVector2d& p2) : P1(p1), P2(p2) {}
};


struct Segment3D {
    Eigen::RowVector3d P1, P2;
    Segment3D(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2) : P1(p1), P2(p2) {}
};

struct Ray3D {
    Eigen::RowVector3d P1, P2;
    Ray3D(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2) : P1(p1), P2(p2) {}
};

struct Line3D {
    Eigen::RowVector3d P1, P2;
    Line3D(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2) : P1(p1), P2(p2) {}
};

#endif
