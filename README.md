# lightGeo

**lightGeo** is a lightweight geometry computation library developed using [Eigen](https://eigen.tuxfamily.org/) for fast 2D/3D computations. It provides basic primitives, advanced bounding‑box algorithms, rasterization tools, path‑planning modules, and a simple camera projection model.

---

## ✨ Features

- **Core Geometry**  
  - 2D/3D points, lines, planes, polygons, and polyhedra  
  - Efficient Eigen‑based linear algebra  

- **Oriented Bounding Box (OBB)**  
  - **PCA Approximation**: Fast, O(n) preprocessing + O(d³) eigen decomposition  
  - **Minimal OBB**: Convex‑hull + rotating‑calipers (2D) and CGAL-based optimization (3D)  

- **Rasterization**  
  - 2D/3D line & surface grid‑fill using Bresenham + hash‑grid acceleration  

- **3D Path Planning**  
  - Two variants of Jump Point Search (JPS) for fast obstacle avoidance in large 3D grids  

- **Camera Projection**  
  - Simple pinhole‑style camera model for transforming 3D world coordinates into image pixels  

## ✨ Features

- 📌 Basic 2D and 3D geometric primitives (points, lines, planes, polyhedra).
- 📐 Advanced 2D/3D OBB (Oriented Bounding Box) algorithms:
  - Fast approximate OBB based on PCA.
  - Minimum OBB using convex hull and rotating calipers.
- 🧱 Rasterization functions for 2D/3D lines, faces, and polyhedra:
  - Efficient hashing and Bresenham-based rasterization.
- 🚀 High-performance 3D JPS (Jump Point Search) path planning algorithms for obstacle avoidance in large 3D grid spaces.
- 🎥 Perspective camera model for projecting 3D points into image coordinates.
---

## 📋 Requirements

- **C++**: 17 or higher
- **CGAL**: 5.6
- **Eigen**: 3.4.90
- **Boost**: 1.87.0
- **Abseil**: 20240116

> ⚠️ Make sure these dependencies are properly configured in your environment. Version mismatches may cause incompatibility.

---

## ⚙️ Installation

```bash
# Clone and add to your project
git clone https://github.com/yourusername/lightGeo.git
cd lightGeo
# Include the headers in your CMakeLists.txt or project settings
