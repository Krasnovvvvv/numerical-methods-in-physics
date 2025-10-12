#ifndef NUMERICAL_METHODS_IN_PHYSICS_ODETASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_ODETASK_H

#pragma once
#include "Labs/Lab4/Tasks/IODETaskStrategy.h"
#include <memory>

class ODETask {
public:
    explicit ODETask(std::unique_ptr<IODETaskStrategy> strategy = nullptr)
        : strategy_(std::move(strategy)) {}

    void setStrategy(std::unique_ptr<IODETaskStrategy> strategy) {
        strategy_ = std::move(strategy);
    }

    void run(const std::function<std::vector<double>(double, const std::vector<double>&)>& rhs,
             const std::vector<double>& y0,
             double t0, double tn,
             double h, double tol,
             unsigned short graphNumber,
             Plotter* plotter,
             std::vector<double> steps = {},
             std::vector<double> tolerances = {},
             std::vector<IODESolver*> solvers = {})
    {
        if (!strategy_)
            throw std::runtime_error("Strategy not set!");
        strategy_->run(rhs, y0, t0, tn, h, tol, graphNumber, plotter, steps, tolerances, solvers);
    }

private:
    std::unique_ptr<IODETaskStrategy> strategy_;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_ODETASK_H