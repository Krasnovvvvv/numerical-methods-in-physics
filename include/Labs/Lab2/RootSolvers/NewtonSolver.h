#ifndef NUMERICAL_METHODS_IN_PHYSICS_NEWTONSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_NEWTONSOLVER_H

#pragma once
#include "Base/IRootSolver.h"

class NewtonSolver : public IRootSolver {
public:
    std::optional<RootSolveResult> solve(
        std::function<double(double)> func,
        std::function<bool(double)> isInDomain,
        double tol,
        double x0,
        double x1 = 1,
        double step = 0.001,
        size_t max_iter = 100
    ) override {
        std::optional<RootSolveResult> result = std::nullopt;
        auto der = [func](double x) {
            const double h = 1e-8;
            return (func(x + h) - func(x - h)) / (2 * h);
        };

        double x_prev = x0, x_next = x0;
        size_t iter = 0;
        std::vector<std::pair<size_t, double>> residuals;
        while (iter < max_iter) {
            double f_val = func(x_prev);
            double f_der = der(x_prev);
            if (std::abs(f_der) < 1e-14) break;
            x_next = x_prev - f_val / f_der;
            residuals.push_back({iter, std::abs(f_val)});
            if (std::abs(x_next - x_prev) < tol) break;
            x_prev = x_next;
            ++iter;
        }
        result = {x_next, iter, std::abs(func(x_next)), residuals};
        return result;
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_NEWTONSOLVER_H