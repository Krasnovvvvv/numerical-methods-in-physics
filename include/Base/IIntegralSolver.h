#ifndef NUMERICAL_METHODS_IN_PHYSICS_IINTEGRALSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_IINTEGRALSOLVER_H

#pragma once
#include <vector>
#include <optional>
#include <functional>

struct IntegrateResult {
    double integral;
    size_t intervals;
    double error;
    std::vector<std::pair<size_t, double>> estimations;
};

class IIntegralSolver {
public:
    virtual ~IIntegralSolver() = default;
    virtual std::optional<IntegrateResult> integrate(
        std::function<double(double)> func,
        double a, double b,
        double tol,
        size_t max_intervals = 1000) = 0;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_IINTEGRALSOLVER_H