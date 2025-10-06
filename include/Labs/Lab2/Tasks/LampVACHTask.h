#ifndef NUMERICAL_METHODS_IN_PHYSICS_LAMPVACHTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_LAMPVACHTASK_H

#pragma once
#include "Base/IRootSolver.h"
#include "Helpers/Plotter.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <numeric>
#include <algorithm>

class LampVACHTask {
public:
    LampVACHTask(IRootSolver& solver, Plotter& plotter)
        : solver(solver), plotter(plotter) {}

    // Оценка кратности на основе последних приближений
    static double estimateMultiplicity(const std::vector<double>& xs) {
        if (xs.size() < 4) return 1.0;
        double d1 = std::abs(xs[xs.size()-3] - xs[xs.size()-4]);
        double d2 = std::abs(xs[xs.size()-2] - xs[xs.size()-3]);
        double d3 = std::abs(xs[xs.size()-1] - xs[xs.size()-2]);
        if (d1 == 0 || d2 == 0) return 1.0;
        double m = std::log(d3/d2) / std::log(d2/d1) + 1.0;
        return m;
    }

    void run(double V_min, double V_max, double dV,
             double R0, double alpha, double T0,
             double sigma, double S, double epsilon,
             double tol = 1e-10, size_t max_iter = 100,
             double Tguess = 500.0,
             const std::string& method_name = "Method")
    {
        std::vector<double> voltages, currents, iterations, func_calls, roots_history;

        for (double V = V_min; V <= V_max + 1e-8; V += dV) {
            auto f = [=](double T) {
                double R = R0 * (1 + alpha * (T - T0));
                double power = epsilon * sigma * S * (std::pow(T, 4) - std::pow(T0, 4));
                return (V * V) / R - power;
            };
            auto isInDomain = [=](double T) { return T > 10 && T < 4000; };

            size_t call_count = 0;
            auto f_counter = [&](double T) {
                ++call_count;
                return f(T);
            };

            auto result = solver.solve(f_counter, isInDomain, tol, Tguess, 0, 0.5, max_iter);

            if (result) {
                double T = result->root;
                double I = V / (R0 * (1 + alpha * (T - T0)));
                voltages.push_back(V);
                currents.push_back(I);
                iterations.push_back(result->iterations);
                func_calls.push_back(static_cast<double>(call_count));
                roots_history.push_back(T);
                std::cout << "V=" << V << ", I=" << I
                          << ", iters=" << result->iterations
                          << ", f-calls=" << call_count << ", T=" << T << std::endl;
            } else {
                voltages.push_back(V);
                currents.push_back(std::nan(""));
                iterations.push_back(std::nan(""));
                func_calls.push_back(std::nan(""));
            }
        }

        double avg_iters = std::accumulate(iterations.begin(), iterations.end(), 0.0) / iterations.size();
        double avg_calls = std::accumulate(func_calls.begin(), func_calls.end(), 0.0) / func_calls.size();
        double multiplicity = estimateMultiplicity(roots_history);

        std::cout << "=== " << method_name << " summary ===\n";
        std::cout << "V_range: [" << V_min << "; " << V_max << "]\n";
        std::cout << "Mean iterations: " << avg_iters
                  << "; Mean function calls: " << avg_calls << std::endl;
        std::cout << "Estimated multiplicity (by last roots): " << multiplicity << std::endl;

        plotter.plot(voltages, currents,
            method_name + " V-I characteristic", "Voltage, V", "Current, A");

    }
private:
    IRootSolver& solver;
    Plotter& plotter;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_LAMPVACHTASK_H