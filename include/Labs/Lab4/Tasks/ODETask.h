#ifndef NUMERICAL_METHODS_IN_PHYSICS_ODETASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_ODETASK_H

#pragma once
#include "Base/IODESolver.h"
#include "Helpers/Plotter.h"
#include "Helpers/Timer.h"
#include <iostream>
#include <optional>

class ODETask {
public:
    explicit ODETask(IODESolver& solver, Plotter* plotter = nullptr, unsigned short graph_num = 1)
        : solver(solver), plotter(plotter), graphNumber(graph_num) {}

void run(
const std::function<std::vector<double>(double, const std::vector<double>&)>& func,
    const std::vector<double>& y0, double t0, double tn, double h, double tol = 1e-8) {
        Timer<std::chrono::microseconds> timer;
        std::optional<ODEResult> result = solver.solve(func, y0, t0, tn, h, tol);
        auto elapsed_us = timer.elapsed();

        std::cout << "\n=== " << solver.name() << " ===\n";

        if (result) {
            std::cout << "Grid step: " << h << "\n"
                      << "Number of steps: " << result->steps << "\n"
                      << "Maximum error per step: "
                      << *std::max_element(result->errorEstimates.begin(),
                          result->errorEstimates.end()) << std::endl;

            if (plotter && !result.y.empty()) {
                std::vector<double> x, v;
                for (auto& state : result.y) {
                    x.emplace_back(state[0]);
                    v.emplace_back(state[1]);
                }
                switch (graphNumber) {
                    case 1:
                        plotter->plot(result.t, x, "x(t)", "t", "x");
                        break;
                    case 2:
                        plotter->plot(result.t, v, "v(t)", "t", "v");
                        break;
                    case 3:
                        plotter->plot(result.t, result.errorEstimates,
                            "Error for " + solver.name(), "t", "error");
                        break;
                    default:
                        break;
                }
            }
       }
    }

private:
    IODESolver& solver;
    Plotter* plotter;
    unsigned short graphNumber;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_ODETASK_H