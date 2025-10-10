#ifndef NUMERICAL_METHODS_IN_PHYSICS_CENTRALRECTSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_CENTRALRECTSOLVER_H

#pragma once
#include "Base/IIntegralSolver.h"
#include <cmath>

class CentralRectSolver : public IIntegralSolver {
public:
    std::string name() const override { return "The middle rectangles"; }
    std::optional<IntegrateResult> integrate(
        std::function<double(double)> func, double a, double b,
        double tol, size_t max_intervals = 1000
    ) override {
        size_t n = 2;
        double prev = 0, curr, h, error = 0;
        std::vector<std::pair<size_t, double>> estim;
        std::vector<double> errors_hist;

        for (; n <= max_intervals; n *= 2) {
            h = (b - a) / n;
            curr = 0;
            for (size_t i = 0; i < n; ++i)
                curr += func(a + h * (i + 0.5));
            curr *= h;
            estim.emplace_back(n, curr);
            if (n > 2) {
                error = std::abs(curr - prev) / 3.0;
                errors_hist.push_back(error);
                if (error < tol)
                    break;
            }
            prev = curr;
        }
        return IntegrateResult{curr, n, error, estim, errors_hist};
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_CENTRALRECTSOLVER_H
