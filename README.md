# LiteGeometry 🌟

**LiteGeometry** is a lightweight geometry computation tool built on [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) for high-performance 2D/3D geometric computations.

---
## ✨ Features

- **🧩 Geometric Primitives**  
  - 2D/3D points, lines, planes, polygons, and polyhedra  
  - Efficient Eigen‑based linear algebra  

- **📦 Oriented Bounding Box (OBB)**  
  - **PCA Approximation**: Fast approximate OBB based on PCA. 
  - **Minimal OBB**: Minimum OBB using convex hull and rotating calipers.

- **🧱 Efficient Rasterization**  
  - 2D/3D line base on Bresenham's algorithm
  - 2D/3D face Rasterization
  - 2D/3D polygon Rasterization
  - polyhedra's face Rasterization 

- **🚀 3D Path Planning**  
  - High-performance 3D JPS (Jump Point Search) path planning algorithms for obstacle avoidance in large 3D map. 

- **📷 Camera Projection**  
  - Perspective camera model for projecting 3D world coordinates into image coordinates  
---

## 🛠 Development Environment
This library is actively developed and tested with:
- **Compiler**: MSVC 2022 (C++17 mode)
- **Eigen**: 3.4.90
- **CGAL**: 5.6
- **Boost**: 1.87.0
- **Abseil**: 20240116

> ⚠️ Earlier versions may work but aren't fully validated.

---

## ⚙️ Installation
# Include the headers in LiteGeometry.h or project settings

## License
LiteGeometry software: [BSD-3-Clause License](LICENSE.txt)
