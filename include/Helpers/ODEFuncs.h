#ifndef NUMERICAL_METHODS_IN_PHYSICS_ODEFUNCS_H
#define NUMERICAL_METHODS_IN_PHYSICS_ODEFUNCS_H

#pragma once
#include <vector>

inline std::vector<double> oscillator(double t, const std::vector<double>& y) {
    constexpr double omega = 10.0;
    constexpr double b = 1.0;
    return {
        y[1],
        - (2 * b / t) * y[1] - omega * omega * y[0]
    };
}

inline std::vector<double> harmonicOscillator(double t, const std::vector<double>& y) {
    constexpr double omega = 1.0;
    return {
        y[1],
        - omega * omega * y[0]
        };
}

#endif //NUMERICAL_METHODS_IN_PHYSICS_ODEFUNCS_H