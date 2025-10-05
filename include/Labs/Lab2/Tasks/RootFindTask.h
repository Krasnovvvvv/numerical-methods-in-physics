#ifndef NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H

#pragma once
#include "Base/IRootSolver.h"
#include "Helpers/Timer.h"
#include <iostream>
#include <iomanip>
#include <functional>
#include <optional>
#include <chrono>

class RootFindTask {
public:
    RootFindTask(IRootSolver& solver) : solver(solver) {}

    void run(
        std::function<double(double)> func,
        std::function<bool(double)> isInDomain,
        double tol,
        double x0,
        double x1 = 1,
        double step = 0.001,
        size_t max_iter = 100
    ) {
        Timer<std::chrono::microseconds> timer;

        std::optional<RootSolveResult> result =
            solver.solve(func, isInDomain, tol, x0, x1, step, max_iter);

        auto elapsed_us = timer.elapsed();
        std::cout << std::setprecision(14);

        if (result) {
            std::cout << "Root: " << result->root << std::endl
                      << "Iterations: " << result->iterations << std::endl
                      << "Residual: " << result->residual << std::endl
                      << "Elapsed: " << elapsed_us << " microseconds" << std::endl;
        } else {
            std::cout << "No root found on specified interval!" << std::endl;
        }
    }

private:
    IRootSolver& solver;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H
