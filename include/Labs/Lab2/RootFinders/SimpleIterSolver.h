#ifndef NUMERICAL_METHODS_IN_PHYSICS_SIMPLEITERSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_SIMPLEITERSOLVER_H

#pragma once
#include "Base/IRootSolver.h"

class SimpleIterSolver : public IRootSolver {
public:
    SimpleIterSolver(double tau) : tau(tau) {}

    RootSolveResult solve(
        std::function<double(double)> func,
        std::function<bool(double)> isInDomain,
        double x0, double tol,
        size_t max_iter = 100
    ) override {
        auto phi = [func, tau](double x) { return x - tau * func(x); };
        double x_prev = x0, x_next = phi(x0);
        size_t iter = 0;
        std::vector<std::pair<size_t, double>> residuals;
        while (std::abs(x_next - x_prev) > tol && iter < max_iter) {
            x_prev = x_next;
            x_next = phi(x_prev);
            ++iter;
            residuals.push_back({iter, std::abs(func(x_next))});
        }
        return {x_next, iter, std::abs(func(x_next)), residuals};
    }

private:
    double tau;
};


#endif //NUMERICAL_METHODS_IN_PHYSICS_SIMPLEITERSOLVER_H