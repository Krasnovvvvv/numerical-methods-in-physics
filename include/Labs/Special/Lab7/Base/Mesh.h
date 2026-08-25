#ifndef NUMERICAL_METHODS_IN_PHYSICS_MESH_H
#define NUMERICAL_METHODS_IN_PHYSICS_MESH_H
#pragma once

#include <stdexcept>

class Mesh {
public:
    int Nx{};
    int Ny{};
    double Lx{};
    double Ly{};
    double dx{};
    double dy{};

    Mesh() = default;

    Mesh(int nx, int ny, double lx, double ly)
        : Nx(nx), Ny(ny), Lx(lx), Ly(ly) {
        if (Nx <= 0 || Ny <= 0) {
            throw std::runtime_error("Mesh: Nx and Ny must be positive.");
        }
        if (Lx <= 0.0 || Ly <= 0.0) {
            throw std::runtime_error("Mesh: Lx and Ly must be positive.");
        }
        dx = Lx / static_cast<double>(Nx);
        dy = Ly / static_cast<double>(Ny);
    }

    int cellCount() const { return Nx * Ny; }
    int id(int i, int j) const { return j * Nx + i; }

    bool inside(int i, int j) const {
        return 0 <= i && i < Nx && 0 <= j && j < Ny;
    }

    bool hasWest(int i) const  { return i > 0; }
    bool hasEast(int i) const  { return i < Nx - 1; }
    bool hasSouth(int j) const { return j > 0; }
    bool hasNorth(int j) const { return j < Ny - 1; }

    double xc(int i) const { return (static_cast<double>(i) + 0.5) * dx; }
    double yc(int j) const { return (static_cast<double>(j) + 0.5) * dy; }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_MESH_H