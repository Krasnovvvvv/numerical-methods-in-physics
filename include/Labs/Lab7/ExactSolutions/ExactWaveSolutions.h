#ifndef EXACT_WAVE_SOLUTIONS_H
#define EXACT_WAVE_SOLUTIONS_H

#define _USE_MATH_DEFINES

#include <cmath>
#include <algorithm>
#include "Labs/Lab7/Tasks/TaskWaveBase.h"

inline double exact_piecewise_velocity(double x, double t,
                                       const TaskWaveBase::PhysParams& phys,
                                       int Nmax = 1000)
{
    double c = phys.c;
    double L = phys.L;
    double x0 = phys.x0;
    double d = phys.delta;
    double v0 = phys.v0;

    double a = std::max(0.0, x0 - d);
    double b = std::min(L, x0 + d);

    if (a >= b) return 0.0;

    double sum = 0.0;
    for (int n = 1; n <= Nmax; ++n) {
        double kn = n * M_PI / L;
        double In = v0 / kn * (std::cos(kn * a) - std::cos(kn * b));
        double Bn = 2.0 / (c * L * kn) * In;
        sum += Bn * std::sin(kn * x) * std::sin(kn * c * t);
    }

    return sum;
}

#endif
