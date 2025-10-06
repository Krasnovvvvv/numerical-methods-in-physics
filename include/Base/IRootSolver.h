#ifndef NUMERICAL_METHODS_IN_PHYSICS_IROOTSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_IROOTSOLVER_H

#pragma once
#include <vector>
#include <functional>
#include <optional>

struct RootSolveResult {
    double root;
    size_t iterations;
    double residual;
    std::vector<std::pair<size_t, double>> residuals;
};

class IRootSolver {
public:
    virtual ~IRootSolver() = default;
    virtual std::optional<RootSolveResult> solve(
        std::function<double(double)> func,
        std::function<bool(double)> isInDomain,
        double tol,
        double x0,
        double x1 = 1, // для дихотомии - правый конец, для односторонних методов игнорируется
        double step = 0.001,
        size_t max_iter = 100
    ) = 0;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_IROOTSOLVER_H