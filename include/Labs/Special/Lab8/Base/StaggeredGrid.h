#ifndef NUMERICAL_METHODS_IN_PHYSICS_STAGGEREDGRID_H
#define NUMERICAL_METHODS_IN_PHYSICS_STAGGEREDGRID_H
#pragma once
#include <cassert>

struct StaggeredGrid {
    int Nx = 0;
    int Ny = 0;
    double Lx = 1.0;
    double Ly = 1.0;

    StaggeredGrid() = default;
    StaggeredGrid(int nx, int ny, double lx, double ly)
        : Nx(nx), Ny(ny), Lx(lx), Ly(ly) {
        assert(Nx > 1);
        assert(Ny > 1);
    }

    double dx() const { return Lx / static_cast<double>(Nx); }
    double dy() const { return Ly / static_cast<double>(Ny); }

    int pNx() const { return Nx; }
    int pNy() const { return Ny; }
    int uNx() const { return Nx + 1; }
    int uNy() const { return Ny; }
    int vNx() const { return Nx; }
    int vNy() const { return Ny + 1; }
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_STAGGEREDGRID_H
