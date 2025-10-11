#ifndef NUMERICAL_METHODS_IN_PHYSICS_INTEGRATETASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_INTEGRATETASK_H

#pragma once
#include "Base/IIntegralSolver.h"
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

class IntegrateTask {
public:
    explicit IntegrateTask(IIntegralSolver& integrator, Plotter* plotter = nullptr)
            : integrator(integrator), plotter(plotter) {}

    void run(
            std::function<double(double)> func,
            double a, double b,
            double tol = 1e-8,
            size_t max_intervals = 4096
    ) {
        Timer<std::chrono::microseconds> timer;
        std::optional<IntegrateResult> result = integrator.integrate(std::move(func), a, b, tol, max_intervals);
        auto elapsed_us = timer.elapsed();

        std::cout << "\n=== " << integrator.name() << " ===\n";
        if (result) {
            std::cout << "Integral: " << std::fixed << std::setprecision(8) << result->integral << "\n"
                      << "Nodes: " << result->intervals << "\n"
                      << "Error estimation: " << std::scientific << result->error << "\n"
                      << "Time: " << elapsed_us << " microseconds" << std::endl;

            if (plotter && result->estimations.size() > 2 && result->errors_hist.size() > 1) {
                // --- Graph 1: ln(error) vs ln(n) ---
                std::vector<double> x_ln, y_lnerr;
                for (size_t i = 1; i < result->estimations.size(); ++i) {
                    x_ln.emplace_back(std::log(static_cast<double>(result->estimations[i].first))); // ln(n)
                    y_lnerr.emplace_back(std::log(std::abs(result->errors_hist[i-1])));            // ln(error)
                }
                plotter->plot(
                    x_ln, y_lnerr,
                    "ln(error) vs ln(n) for " + integrator.name(),
                    "ln(n)",
                    "ln(error)",
                    false);

                // --- Graph 2: gamma_k vs ln(n) ---
                std::vector<double> x, y;
                auto safe_log = [](double v) { return v <= 0 ? NAN : std::log(v); };
                for (size_t i = 1; i + 1 < result->estimations.size(); ++i) {
                    double ln_nk1 = safe_log(result->estimations[i + 1].first);
                    double ln_nk  = safe_log(result->estimations[i].first);
                    double ln_Rk1 = safe_log(std::abs(result->errors_hist[i]));
                    double ln_Rk  = safe_log(std::abs(result->errors_hist[i - 1]));
                    if (std::isnan(ln_nk1) || std::isnan(ln_nk) || std::isnan(ln_Rk) || std::isnan(ln_Rk1))
                        continue;
                    double denom = ln_nk1 - ln_nk;
                    if (std::abs(denom) < 1e-8) continue;
                    double gamma_k = (ln_Rk - ln_Rk1) / denom;
                    x.emplace_back((ln_nk + ln_nk1) / 2.0);
                    y.emplace_back(gamma_k);
                }

                plotter->plot(
                        x, y,
                        "gamma_k vs ln(n) for " + integrator.name(),
                        "ln(n)",
                        "gamma_k",
                        false);
            }
        } else {
            std::cout << "Integral is not calculated (error)\n";
        }
    }

private:
    IIntegralSolver& integrator;
    Plotter* plotter;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_INTEGRATETASK_H
