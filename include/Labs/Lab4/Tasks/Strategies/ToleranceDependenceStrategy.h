#ifndef NUMERICAL_METHODS_IN_PHYSICS_TOLERANCEDEPENDENCESTRATEGY_H
#define NUMERICAL_METHODS_IN_PHYSICS_TOLERANCEDEPENDENCESTRATEGY_H

#pragma once
#include "Base/IODESolver.h"
#include "Labs/Lab4/Tasks/IODETaskStrategy.h"

class ToleranceDependenceStrategy : public IODETaskStrategy {
public:
    ToleranceDependenceStrategy(IODESolver* solver) : solver_(solver) {}

    void run(const std::function<std::vector<double>(double, const std::vector<double>&)>& rhs,
             const std::vector<double>& y0,
             double t0, double tn,
             double h, double /*tol*/,
             unsigned short graphNumber,
             Plotter* plotter,
             std::vector<double> /*steps*/,
             std::vector<double> tolerances,
             std::vector<IODESolver*> /*solvers*/) override
    {
        if (!solver_ || !plotter) throw std::runtime_error("Solver or Plotter not set!");
        if (tolerances.empty()) throw std::invalid_argument("Tolerances vector is empty!");

        std::vector<std::vector<double>> xs, ys;
        std::vector<std::string> labels;

        for (double tol : tolerances) {
            auto res = solver_->solve(rhs, y0, t0, tn, h, tol);

            std::vector<double> x, v;
            for (const auto& state : res.y) {
                x.push_back(state[0]);
                v.push_back(state[1]);
            }

            switch (graphNumber) {
                case 1: xs.push_back(res.t); ys.push_back(x); break;
                case 2: xs.push_back(res.t); ys.push_back(v); break;
                case 3: xs.push_back(res.t); ys.push_back(res.errorEstimates); break;
                case 4: xs.push_back(x); ys.push_back(v); break;
                default: break;
            }
            labels.push_back("tol = " + std::to_string(tol));
        }
        plotter->plot(xs, ys, labels,
            (graphNumber == 4 ? "x" : "t"),
            (graphNumber == 4 ? "v" : (graphNumber == 1 ? "x" : graphNumber == 2 ? "v" : "error")));
    }

private:
    IODESolver* solver_;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_TOLERANCEDEPENDENCESTRATEGY_H