#ifndef NUMERICAL_METHODS_IN_PHYSICS_SCALARFIELD_H
#define NUMERICAL_METHODS_IN_PHYSICS_SCALARFIELD_H
#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include "Mesh.h"

class ScalarField {
public:
    Mesh mesh;
    std::vector<double> data;

    ScalarField() = default;

    explicit ScalarField(const Mesh& m, double initialValue = 0.0)
        : mesh(m), data(static_cast<size_t>(m.cellCount()), initialValue) {}

    double& operator()(int i, int j) {
        return data.at(static_cast<size_t>(mesh.id(i, j)));
    }

    const double& operator()(int i, int j) const {
        return data.at(static_cast<size_t>(mesh.id(i, j)));
    }

    int size() const {
        return static_cast<int>(data.size());
    }

    void fill(double value) {
        std::fill(data.begin(), data.end(), value);
    }

    double maxAbsDiff(const ScalarField& other) const {
        if (mesh.Nx != other.mesh.Nx || mesh.Ny != other.mesh.Ny ||
            mesh.Lx != other.mesh.Lx || mesh.Ly != other.mesh.Ly) {
            throw std::runtime_error("ScalarField::maxAbsDiff: incompatible meshes.");
            }

        double mx = 0.0;
        for (size_t k = 0; k < data.size(); ++k) {
            mx = std::max(mx, std::abs(data[k] - other.data[k]));
        }
        return mx;
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_SCALARFIELD_H