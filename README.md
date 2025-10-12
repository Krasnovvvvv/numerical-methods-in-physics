# ⚡ Numerical Methods in Physics

![GitHub last commit](https://img.shields.io/github/last-commit/Krasnovvvvv/numerical-methods-in-physics)
![GitHub top language](https://img.shields.io/github/languages/top/Krasnovvvvv/numerical-methods-in-physics)
![MIT License](https://img.shields.io/github/license/Krasnovvvvv/numerical-methods-in-physics)
![Issues](https://img.shields.io/github/issues/Krasnovvvvv/numerical-methods-in-physics)
![Stars](https://img.shields.io/github/stars/Krasnovvvvv/numerical-methods-in-physics)

![Ubuntu](https://img.shields.io/badge/Ubuntu-20.04+-orange?logo=ubuntu)
![Windows](https://img.shields.io/badge/Windows-10+-blue?logo=windows)
![CI Ubuntu](https://img.shields.io/github/actions/workflow/status/Krasnovvvvv/numerical-methods-in-physics/ci.yml?branch=main&label=Ubuntu&logo=ubuntu)
![CI Windows](https://img.shields.io/github/actions/workflow/status/Krasnovvvvv/numerical-methods-in-physics/ci.yml?branch=main&label=Windows&logo=windows)

This project contains C++ implementations of a variety of classic and modern numerical techniques, applied to physical problems encountered in laboratory courses and practical assignments. The codebase is designed with clarity, scalability, and reproducibility in mind

---

## ℹ️ About

These laboratory works showcase algorithms for solving core computational physics problems. The implementations are modular and leverage clean C++ patterns to facilitate study, reuse, and extension for new problems

---

## 🔬 Laboratory Topics

### 🧮 Systems of Linear Algebraic Equations (SLAE)

Demonstrates solving tridiagonal and arbitrary systems using direct and iterative methods (e.g., Thomas algorithm, Jacobi, Gauss-Seidel), with application to real modeling tasks

### 🌱 Root Finding Algorithms 

Contains implementations of the bisection (dichotomy), Newton, and simple iteration methods for solving non-linear equations. Demonstrates convergence efficiency and robustness on real physical functions

### 📐 Numerical Integration 

Features various quadrature methods (rectangular, trapezoidal, Simpson's), with error analysis and comparisons on benchmark test cases

### 📊 Ordinary Differential Equations (ODEs) 

Includes explicit and implicit schemes for ODEs, such as Euler, improved Euler and Runge-Kutta methods, tested on classical physical systems (oscillator, decay, etc.)

---

## 🚀 Getting Started 

### 🛠️ Prerequisites

- C++17 compatible compiler (GCC, Clang, MSVC)
- CMake 3.14 or newer
- [matplot++](https://alandefreitas.github.io/matplotplusplus/) development files
- [Eigen](https://eigen.tuxfamily.org/) library (header-only, managed via CMake or package manager)

### ⚡ Build

To build the code, standard CMake workflows are recommended:

```bash
git clone https://github.com/Krasnovvvvv/numerical-methods-in-physics.git
cd numerical-methods-in-physics
mkdir build && cd build
cmake ..
make
./lab_executable
```

---

## 📄 Reports

All completed laboratory reports (including results, and code descriptions) can be found here:

|        📝 Report                         |           Link                |
|:----------------------------------------:|:-----------------------------:|
| Systems of Linear Algebraic Equations 📄 | [Read](reports/Lab1.md)       |
| Root Finding Algorithms 📄               | [Read](reports/Lab2.md)       |

---

## 📦 Dependencies

- [matplot++](https://alandefreitas.github.io/matplotplusplus/) — a high-quality C++ plotting library for scientific visualization
- [Eigen](https://eigen.tuxfamily.org/) — a fast, versatile C++ library for linear algebra and matrix operations

> These dependencies are automatically handled via CMake (assuming installed on your system or via package managers like vcpkg)

---

## 🤝 Contribution

Contributions, bug reports, and feature requests are welcome!  
Please open issues or pull requests on the GitHub repository

---

## 📝 License

This project is licensed under the MIT License — see [![License](https://img.shields.io/github/license/Krasnovvvvv/numerical-methods-in-physics)](LICENSE)






