#ifndef NUMERICAL_METHODS_IN_PHYSICS_IODETASKSTRATEGY_H
#define NUMERICAL_METHODS_IN_PHYSICS_IODETASKSTRATEGY_H

#pragma once
#include <functional>
#include "Base/IODESolver.h"
#include "Helpers/Plotter.h"

class IODETaskStrategy {
public:
    virtual ~IODETaskStrategy() = default;
    virtual void run(
        const std::function<std::vector<double>(double, const std::vector<double>&)>& rhs,
        const std::vector<double>& y0,
        double t0, double tn,
        double h, double tol,
        unsigned short graphNumber,
        Plotter* plotter,
        std::vector<double> steps = {},
        std::vector<double> tolerances = {},
        std::vector<IODESolver*> solvers = {}
    ) = 0;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_IODETASKSTRATEGY_H