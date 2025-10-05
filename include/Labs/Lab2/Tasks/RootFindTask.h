#ifndef NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H

#pragma once
#include "Base/IRootFinder.h"
#include "Helpers/Timer.h"
#include <iostream>
#include <iomanip>
#include <functional>
#include <chrono>

class RootFindTask {
public:
    RootFindTask(IRootFinder& finder) : finder(finder) {}

    void run(
        std::function<double(double)> func,
        std::function<bool(double)> isInDomain,
        double a, double b,
        double step, double tol
    ) {
        Timer<std::chrono::microseconds> timer;
        auto result = finder.findSignChange(func, isInDomain, a, b, step, tol);
        auto elapsed_us = timer.elapsed();

        if (result) {
            std::cout << "Interval of sign change: ["
                      << result->left << ", " << result->right << "]\n";
            std::cout << "f(left) = " << func(result->left)
                      << ", f(right) = " << func(result->right) << std::endl;
            std::cout << "Elapsed: " << (elapsed_us)
                      << " microseconds (" << elapsed_us << " us)\n";
        } else {
            std::cout << "No sign change found in interval ["
                      << a << ", " << b << "]\n";
        }
    }

private:
    IRootFinder& finder;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H
