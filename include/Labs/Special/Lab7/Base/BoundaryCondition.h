#ifndef NUMERICAL_METHODS_IN_PHYSICS_BOUNDARYCONDITION_H
#define NUMERICAL_METHODS_IN_PHYSICS_BOUNDARYCONDITION_H
#pragma once

#include <functional>
#include <utility>

enum class BoundaryType {
    Dirichlet,
    OutletZeroGradient
};

struct BoundaryCondition {
    BoundaryType type = BoundaryType::Dirichlet;
    std::function<double(double)> value = [](double) { return 0.0; };

    BoundaryCondition() = default;

    BoundaryCondition(BoundaryType t, std::function<double(double)> v)
        : type(t), value(std::move(v)) {}

    static BoundaryCondition dirichlet(double c) {
        return BoundaryCondition(
            BoundaryType::Dirichlet,
            [c](double) { return c; }
        );
    }

    static BoundaryCondition dirichlet(const std::function<double(double)>& f) {
        return BoundaryCondition(BoundaryType::Dirichlet, f);
    }

    static BoundaryCondition outlet() {
        return BoundaryCondition(
            BoundaryType::OutletZeroGradient,
            [](double) { return 0.0; }
        );
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_BOUNDARYCONDITION_H