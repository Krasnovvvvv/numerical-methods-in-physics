#ifndef NUMERICAL_METHODS_IN_PHYSICS_CYLINDERHEATPROBLEM_H
#define NUMERICAL_METHODS_IN_PHYSICS_CYLINDERHEATPROBLEM_H

#pragma once
#include <functional>
#include <cstddef>

struct CylinderHeatProblem {
    size_t nr = 50;    // число внутренних узлов по r
    size_t nz = 50;    // число внутренних узлов по z
    double H = 1.0;    // H = l / a
    double q = 0.0;    // dT/dz |_{z=0} = -q
    double theta = 0.0;// T |_{z=H} = theta

    std::function<double(double, double)> rhs =
        [](double, double) { return 0.0; };
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_CYLINDERHEATPROBLEM_H
