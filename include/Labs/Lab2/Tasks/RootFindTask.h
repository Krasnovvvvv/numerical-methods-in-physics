#ifndef NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H

#pragma once
#include "Base/IRootFinder.h"
#include "Helpers/Timer.h"
#include <iostream>
#include <iomanip>
#include <functional>
#include <chrono>
#include <iomanip>

class RootFindTask {
public:
    RootFindTask(IRootFinder& finder) : finder(finder) {}

    void run(
        std::function<double(double)> func,
        std::function<bool(double)> isInDomain,
        double a, double b,
        double step, double tol, size_t max_iter = 100
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
                      << " microseconds\n";
            // Уточняем корень методом дихотомии на найденном интервале:
            double root = refineRootBisection(func, result->left, result->right, tol);
            std::cout << "Refined root: x = " <<std::setprecision(15)<< root
                      << ", f(x) = " << func(root) << std::endl;
        } else {
            std::cout << "No sign change found in interval ["
                      << a << ", " << b << "]\n";
        }

    }

private:
    IRootFinder& finder;

    double refineRootBisection(
    std::function<double(double)> f,
    double a, double b,
    double tol = 1e-12,
    size_t max_iter = 100
) {
        double fa = f(a), fb = f(b);
        if (fa * fb > 0)
            throw std::runtime_error("No sign change on interval for bisection!");

        for (size_t iter = 0; iter < max_iter && (b - a) > tol; ++iter) {
            double c = 0.5 * (a + b);
            double fc = f(c);
            if (std::abs(fc) < tol)
                return c;
            if (fa * fc < 0) {
                b = c;
                fb = fc;
            } else {
                a = c;
                fa = fc;
            }
        }
        return 0.5 * (a + b);
    }

};

#endif //NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H
