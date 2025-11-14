#ifndef NUMERICAL_METHODS_IN_PHYSICS_ODETASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_ODETASK_H

#pragma once
#include "Base/IODESolver.h"
#include "Helpers/Timer.h"
#include "Helpers/Plotter.h"
#include <iostream>
#include <iomanip>
#include <functional>
#include <optional>
#include <chrono>
#include <utility>
#include <cmath>
#include <vector>

class ODETask {
public:
    explicit ODETask(IODESolver& solver, Plotter* plotter = nullptr, unsigned short graph_num = 1)
            : solver(solver), plotter(plotter), graphNumber(graph_num) {}

    void run(
        std::function<std::vector<double>(double, const std::vector<double>&)> func,
        std::vector<double> y0,
        double t0, double tn, double h,
        std::function<double(double)> exact_sol = nullptr
    ) {
        Timer<std::chrono::microseconds> timer;
        auto result = solver.solve(func, y0, t0, tn, h);
        auto elapsed_us = timer.elapsed();

        std::cout << "\n=== " << solver.name() << " ===\n";
        std::cout << "Steps: " << result.steps << "\n";
        std::cout << "Computation time: " << elapsed_us << " microseconds\n";
        std::cout << std::fixed << std::setprecision(7);
        std::cout << " t      | Num. U    | Error\n";
        for (size_t i = 0; i < result.t.size(); ++i) {
            double exact = exact_sol ? exact_sol(result.t[i]) : NAN;
            double error = (exact_sol && i < result.y.size()) ? std::abs(result.y[i][0] - exact) : NAN;
            std::cout << std::setw(7) << result.t[i] << " | "
                      << std::setw(9) << result.y[i][0] << " | "
                      << std::setw(9) << error << "\n";
        }

        // ------ Graphical output ------
        if (plotter) {
            std::vector<double> x, y_numeric, y_exact, y_err;
            for (size_t i = 0; i < result.t.size(); ++i) {
                x.push_back(result.t[i]);
                y_numeric.push_back(result.y[i][0]);
                if (exact_sol) y_exact.push_back(exact_sol(result.t[i]));
            }
            if (exact_sol) {
                for (size_t i = 0; i < y_numeric.size(); ++i)
                    y_err.push_back(std::abs(y_numeric[i] - y_exact[i]));
            }

            if (graphNumber == 1) {
                // Numerical and exact solution on one graph
                std::vector<std::vector<double>> all_x, all_y;
                std::vector<std::string> labels;
                all_x.push_back(x);
                all_y.push_back(y_numeric);
                labels.push_back("Numerical solution: " + solver.name());
                if (exact_sol) {
                    all_x.push_back(x);
                    all_y.push_back(y_exact);
                    labels.push_back("Exact solution");
                }
                plotter->plot(all_x, all_y, labels, "t", "U");
            }
            if (graphNumber == 2 && exact_sol) {
                // Error graph
                plotter->plot(x, y_err, "Error for: " + solver.name(), "t", "Error");
            }
        }
    }

private:
    IODESolver& solver;
    Plotter* plotter;
    unsigned short graphNumber;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_ODETASK_H
