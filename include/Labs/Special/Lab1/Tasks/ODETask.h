#ifndef NUMERICAL_METHODS_IN_PHYSICS_ODETASK_H
#define NUMERIСAL_METHODS_IN_PHYSICS_ODETASK_H

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
#include <string>

namespace special {
    class ODETask {
    public:
        explicit ODETask(
            IODESolver& solver,
            Plotter* plotter = nullptr,
            unsigned short graph_num = 1,
            std::vector<std::string> component_names = {},
            std::string time_name = "t"
        )
        : solver(solver), plotter(plotter), graphNumber(graph_num),
          componentNames(std::move(component_names)), timeName(std::move(time_name)) {}

        void run(
            std::function<std::vector<double>(double, const std::vector<double>&)> func,
            std::vector<double> y0,
            double t0, double tn, double h,
            std::vector<std::function<double(double)>> exact_solutions = {}
        ) {
            Timer<std::chrono::microseconds> timer;
            auto result = solver.solve(func, y0, t0, tn, h);
            auto elapsed_us = timer.elapsed();

            const size_t dim = result.y[0].size();
            std::cout << "\n=== " << solver.name() << " ===\n";
            std::cout << "Steps: " << result.steps << "\n";
            std::cout << "Computation time: " << elapsed_us << " microseconds\n";
            std::cout << std::fixed << std::setprecision(7);

            // --- Euclidean error for each time point ---
            std::vector<double> euclid_errors(result.t.size(), NAN);
            if (!exact_solutions.empty()) {
                for (size_t i = 0; i < result.t.size(); ++i) {
                    double sum_sq = 0.0;
                    for (size_t comp = 0; comp < dim && comp < exact_solutions.size(); ++comp) {
                        double val_num = result.y[i][comp];
                        double val_ex  = exact_solutions[comp](result.t[i]);
                        double delta = val_num - val_ex;
                        sum_sq += delta * delta;
                    }
                    euclid_errors[i] = std::sqrt(sum_sq);
                }
            }

            for (size_t comp = 0; comp < dim; ++comp) {
                std::string comp_name = (comp < componentNames.size()) ?
                componentNames[comp] : ("U" + std::to_string(comp+1));
                std::cout << "\n--- " << comp_name << " ---\n";
                std::cout << " t      | Num. U    | Error\n";
                for (size_t i = 0; i < result.t.size(); ++i) {
                    double val_num = result.y[i][comp];
                    double val_ex = (!exact_solutions.empty() && comp < exact_solutions.size()) ?
                        exact_solutions[comp](result.t[i]) : NAN;
                    double err = (!exact_solutions.empty() && comp < exact_solutions.size()) ?
                        std::abs(val_num - val_ex) : NAN;
                    std::cout << std::setw(7) << result.t[i] << " | "
                              << std::setw(9) << val_num << " | "
                              << std::setw(9) << err << "\n";
                }
            }

            if (plotter) {
                std::vector<double> x(result.t.begin(), result.t.end());

                // Case 1: all solutions and accurate on one graph
                if (graphNumber == 1) {
                    std::vector<std::vector<double>> all_ys;
                    std::vector<std::string> labels;
                    // Numeric solutions
                    for (size_t comp = 0; comp < dim; ++comp) {
                        std::vector<double> y_num;
                        for (size_t i = 0; i < result.t.size(); ++i)
                            y_num.push_back(result.y[i][comp]);
                        all_ys.push_back(y_num);
                        std::string comp_name = (comp < componentNames.size())
                        ? componentNames[comp] : ("U" + std::to_string(comp+1));
                        labels.push_back(comp_name + " numeric by " + solver.name());
                    }
                    // Exact solutions, if provided
                    if (!exact_solutions.empty()) {
                        for (size_t comp = 0; comp < dim; ++comp) {
                            if (comp < exact_solutions.size()) {
                                std::vector<double> y_ex;
                                for (auto t : x)
                                    y_ex.push_back(exact_solutions[comp](t));
                                all_ys.push_back(y_ex);
                                std::string comp_name = (comp < componentNames.size())
                                ? componentNames[comp] : ("U" + std::to_string(comp+1));
                                labels.push_back(comp_name + " exact");
                            }
                        }
                    }
                    std::vector<std::vector<double>> xs(all_ys.size(), x);
                    plotter->plot(xs, all_ys, labels, timeName, "Components");
                }

                // Case 2: graph of the Euclidean error rate
                if (graphNumber == 2 && !exact_solutions.empty()) {
                    std::vector<std::vector<double>> all_errors{euclid_errors};
                    std::vector<std::string> labels{"Euclidean error by " + solver.name()};
                    std::vector<std::vector<double>> xs{ x };
                    plotter->plot(xs, all_errors, labels, timeName, "Euclidean Error");
                }
            }
        }

    private:
        IODESolver& solver;
        Plotter* plotter;
        unsigned short graphNumber;
        std::vector<std::string> componentNames;
        std::string timeName;
    };
} // namespace special

#endif // NUMERICAL_METHODS_IN_PHYSICS_ODETASK_H
