#ifndef NUMERICAL_METHODS_IN_PHYSICS_CLASSICSTRATEGY_H
#define NUMERICAL_METHODS_IN_PHYSICS_CLASSICSTRATEGY_H

#pragma once
#include "Base/IODESolver.h"
#include "Labs/Lab4/Tasks/IODETaskStrategy.h"
#include "Helpers/Timer.h"

class ClassicStrategy : public IODETaskStrategy {
public:
    using solver = IODESolver;
    ClassicStrategy(solver *solver) : solver_(solver) {}

    void run(const std::function<std::vector<double>(double, const std::vector<double>&)>& rhs,
             const std::vector<double>& y0,
             double t0, double tn,
             double h, double tol,
             unsigned short graphNumber,
             Plotter* plotter,
             std::vector<double> /*steps*/,
             std::vector<double> /*tolerances*/,
             std::vector<solver *> /*solvers*/) override
    {
        if (!solver_ || !plotter)
            throw std::runtime_error("Solver or Plotter not set!");

        Timer timer;
        std::optional<ODEResult> result = solver_->solve(rhs, y0, t0, tn, h, tol);
        auto elapsed_us = timer.elapsed();

        if (!result)
            throw std::runtime_error("Solver returned no result!");

        std::vector<double> x, v;
        for (const auto& state : result->y) {
            x.push_back(state[0]);
            v.push_back(state[1]);
        }

        std::cout << "\n=== " << solver_->name() << " ===\n"
                  << "Grid step: " << h << "\n"
                  << "Number of steps: " << result->steps << "\n"
                  << "Maximum error per step: "
                  << *std::max_element(result->errorEstimates.begin(),
                          result->errorEstimates.end()) << std::endl
                  << "Time: " << elapsed_us << " ms"<<std::endl;

        switch (graphNumber) {
            case 1:
                plotter->plot(result->t, x, "x(t) for " +
                    solver_->name(), "t", "x");
                break;
            case 2:
                plotter->plot(result->t, v, "v(t) for " +
                    solver_->name(), "t", "v");
                break;
            case 3:
                plotter->plot(result->t, result->errorEstimates,
                    "Error for " + solver_->name(), "t", "error");
                break;
            case 4:
                plotter->plot(x, v, "Phase trajectory for " +
                    solver_->name(), "x", "v");
                break;
            default:
                break;
        }
    }

private:
    solver *solver_;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_CLASSICSTRATEGY_H