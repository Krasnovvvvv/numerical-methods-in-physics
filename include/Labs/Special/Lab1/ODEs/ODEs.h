#ifndef NUMERICAL_METHODS_IN_PHYSICS_ODES_H
#define NUMERICAL_METHODS_IN_PHYSICS_ODES_H

#pragma once
#include <vector>
#include <cmath>

namespace rigidSystem {
    inline std::vector<double> rigidSystem(double t, const std::vector<double>& Y) {
        double y = Y[0], z = Y[1];
        return std::vector<double>{
            998*y + 1998*z,
            -999*y - 1999*z
        };
    }

    inline double rigidSystemYSolution(double t) {
        return 4 * std::exp(-t) - 3 * std::exp(-1000 * t);
    }

    inline double rigidSystemZSolution(double t) {
        return -2 * std::exp(-t) + 3 * std::exp(-1000 * t);
    }
} // namespace rigidSystem

// U' = 1 + 0.5*U^2
inline std::vector<double> ode_rhs(double t, const std::vector<double>& U) {
    return std::vector<double>{ 1.0 + 0.5 * U[0] * U[0] };
}
#endif //NUMERICAL_METHODS_IN_PHYSICS_ODES_H