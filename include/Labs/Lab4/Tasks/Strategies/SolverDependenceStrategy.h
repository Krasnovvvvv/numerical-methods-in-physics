#ifndef NUMERICAL_METHODS_IN_PHYSICS_SOLVERDEPENDENCESTRATEGY_H
#define NUMERICAL_METHODS_IN_PHYSICS_SOLVERDEPENDENCESTRATEGY_H

#pragma once
#include "Base/IODESolver.h"
#include "Labs/Lab4/Tasks/IODETaskStrategy.h"

class SolverDependenceStrategy : public IODETaskStrategy {
public:
    void run(const std::function<std::vector<double>(double, const std::vector<double>&)>& rhs,
             const std::vector<double>& y0,
             double t0, double tn,
             double h, double tol,
             unsigned short graphNumber,
             Plotter* plotter,
             std::vector<double> /*steps*/,
             std::vector<double> /*tolerances*/,
             std::vector<IODESolver*> solvers) override
    {
        if (solvers.empty()) throw std::invalid_argument("Solvers vector is empty!");
        if (!plotter) throw std::runtime_error("Plotter not set!");

        std::vector<std::vector<double>> xs, ys;
        std::vector<std::string> labels;

        for (auto* solver : solvers) {
            auto res = solver->solve(rhs, y0, t0, tn, h, tol);

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
            labels.push_back(solver->name());
        }
        plotter->plot(xs, ys, labels,
            (graphNumber == 4 ? "x" : "t"),
            (graphNumber == 4 ? "v" : (graphNumber == 1 ? "x" : graphNumber == 2 ? "v" : "error")));
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_SOLVERDEPENDENCESTRATEGY_H