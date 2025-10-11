#ifndef NUMERICAL_METHODS_IN_PHYSICS_IODESOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_IODESOLVER_H

#pragma once
#include <vector>
#include <functional>
#include <string>

struct ODEResult {
    std::vector<double> t;
    std::vector<std::vector<double>> y;
};

class IODESolver {
public:
    virtual ~IODESolver() = default;
    virtual ODEResult solve(
        const std::function<std::vector<double>(double, const std::vector<double>&)>& func,
        std::vector<double> y0, double t0, double tn, double h) = 0;

    virtual std::string name() const = 0;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_IODESOLVER_H