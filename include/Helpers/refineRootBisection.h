#ifndef NUMERICAL_METHODS_IN_PHYSICS_REFINEROOTBISECTION_H
#define NUMERICAL_METHODS_IN_PHYSICS_REFINEROOTBISECTION_H

#pragma once
#include <cmath>
#include <functional>
#include <stdexcept>

inline double refineRootBisection(
    std::function<double(double)> f,
    double a, double b,
    double tol = 1e-12,
    size_t max_iter = 100
) {
    double fa = f(a), fb = f(b);
    if (fa * fb > 0)
        throw std::runtime_error("No sign change on interval for bisection!");

    for (size_t iter = 0; iter < max_iter && (b - a) > tol; ++iter) {
        double c = 0.5 * (a + b);
        double fc = f(c);
        if (std::abs(fc) < tol)
            return c;
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return 0.5 * (a + b);
}

#endif //NUMERICAL_METHODS_IN_PHYSICS_REFINEROOTBISECTION_H