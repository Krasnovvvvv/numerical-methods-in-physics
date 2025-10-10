#ifndef NUMERICAL_METHODS_IN_PHYSICS_TRAPEZOIDSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_TRAPEZOIDSOLVER_H

#pragma once
#include "Base/IIntegralSolver.h"
#include <cmath>

class TrapezoidSolver : public IIntegralSolver {
public:
    std::optional<IntegrateResult> integrate(
            std::function<double(double)> func,
            double a, double b,
            double tol = 1e-8,
            size_t max_intervals = 1000) override {
        size_t n = 2;
        double prev, curr, h, error = 0;
        std::vector<std::pair<size_t, double>> estim;
        for (; n <= max_intervals; n *= 2) {
            h = (b - a) / n;
            curr = 0.5 * (func(a) + func(b));
            for (size_t i = 1; i < n; ++i) {
                curr += func(a + i * h);
            }
            curr *= h;
            estim.emplace_back(n, curr);
            if (n > 2) {
                error = fabs(curr - prev) / 3.0;
                if (error < tol) break;
            }
            prev = curr;
        }
        return IntegrateResult{curr, n, error, estim};
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_TRAPEZOIDSOLVER_H
