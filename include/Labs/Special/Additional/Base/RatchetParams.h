#ifndef NUMERICAL_METHODS_IN_PHYSICS_RATCHETPARAMS_H
#define NUMERICAL_METHODS_IN_PHYSICS_RATCHETPARAMS_H

#pragma once
#include <cmath>

struct RatchetParams {
    double V1;      // V1/kT
    double V2;      // V2/kT
    double alpha;   // f-
    double epsilon; // tau_c/tau_D
    double a;       // dichotomic noise amplitude
    double dt;

    RatchetParams()
        : V1(0.0),
          V2(0.0),
          alpha(0.0),
          epsilon(0.0),
          a(0.0),
          dt(0.0) {}

    double dVdx(const double x) const {
        static constexpr double TWO_PI  = 2.0 * M_PI;
        static constexpr double FOUR_PI = 4.0 * M_PI;
        return TWO_PI * V1 * std::cos(TWO_PI * x)
             + FOUR_PI * V2 * std::cos(FOUR_PI * x);
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_RATCHETPARAMS_H
