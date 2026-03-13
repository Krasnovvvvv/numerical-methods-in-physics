#ifndef NUMERICAL_METHODS_IN_PHYSICS_RATCHETPARAMS_H
#define NUMERICAL_METHODS_IN_PHYSICS_RATCHETPARAMS_H

#pragma once

#include <cmath>

struct RatchetParams {
    double V1 = 0.0;      // V1/kT
    double V2 = 0.0;      // V2/kT
    double alpha = 0.0;   // f-
    double epsilon = 0.0; // tau_c/tau_D
    double a = 0.0;       // dichotomic noise amplitude
    double dt = 0.0;

    RatchetParams() = default;

    [[nodiscard]] double dVdx(const double x) const noexcept {
        static constexpr double TWO_PI = 2.0 * M_PI;
        static constexpr double FOUR_PI = 4.0 * M_PI;
        return TWO_PI * V1 * std::cos(TWO_PI * x)
             + FOUR_PI * V2 * std::cos(FOUR_PI * x);
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_RATCHETPARAMS_H
