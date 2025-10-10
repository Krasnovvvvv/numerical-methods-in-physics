#ifndef NUMERICAL_METHODS_IN_PHYSICS_INTEGRATETASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_INTEGRATETASK_H

#pragma once
#include "Base/IIntegralSolver.h"
#include "Helpers/Timer.h"
#include <iostream>
#include <iomanip>
#include <functional>
#include <optional>
#include <chrono>

class IntegrateTask {
public:
    IntegrateTask(IIntegralSolver& integrator, const std::string& methodName)
            : integrator(integrator), methodName(methodName) {}

    void run(
            std::function<double(double)> func,
            double a, double b,
            double tol = 1e-8,
            size_t max_intervals = 1000
    ) {
        Timer<std::chrono::microseconds> timer;
        std::optional<IntegrateResult> result = integrator.integrate(func, a, b, tol, max_intervals);
        auto elapsed_us = timer.elapsed();

        std::cout << "\n=== " << methodName << " ===\n";
        if (result) {
            std::cout << "Integral: " << std::fixed << std::setprecision(8) << result->integral << "\n"
                      << "Intervals: " << result->intervals << "\n"
                      << "Error estimation: " << std::scientific << result->error << "\n"
                      << "Time: " << elapsed_us << " microseconds" << std::endl;
        } else {
            std::cout << "Integral is not calculated (error)\n";
        }
    }

private:
    IIntegralSolver& integrator;
    std::string methodName;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_INTEGRATETASK_H
