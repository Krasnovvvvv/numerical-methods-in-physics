#ifndef NUMERICAL_METHODS_IN_PHYSICS_IROOTSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_IROOTSOLVER_H

#pragma once
#include <vector>
#include <functional>

struct RootSolveResult {
    double root;
    size_t iterations;
    double residual;
    std::vector<std::pair<size_t, double>> residuals;
};

class IRootSolver {
public:
    virtual ~IRootSolver() = default;
    virtual RootSolveResult solve(
        std::function<double(double)> func,
        std::function<bool(double)> isInDomain,
        double x0, double tol,
        size_t max_iter = 100
    ) = 0;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_IROOTSOLVER_H