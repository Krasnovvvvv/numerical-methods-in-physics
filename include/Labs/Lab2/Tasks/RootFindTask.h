#ifndef NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H

#pragma once
#include "Base/IRootFinder.h"
#include "Helpers/Timer.h"
#include "Helpers/Plotter.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>
#include <limits>
#include <cmath>

class RootFinderTask {
public:
    RootFinderTask(IRootFinder& finder, Plotter* plot = nullptr)
        : finder(finder), plotter(plot) {}

    void run(
        std::function<double(double)> func,
        std::function<bool(double)> isInDomain,
        double a, double b,
        double step, double tol
    ) {
        Timer<> timer;
        auto result = finder.findSignChange(func, isInDomain, a, b, step, tol);
        auto elapsed_us = timer.elapsed();

        // Собираем функцию на сетке
        std::vector<double> x, y;
        for (double val = a; val <= b; val += step) {
            if (isInDomain(val)) {
                x.push_back(val);
                y.push_back(func(val));
            }
        }

        // Функция для expected: выдаёт только значения на найденном интервале, остальные -- nan
        std::optional<Plotter::ExpectedCurveFunc> mark_segment = std::nullopt;
        if (result) {
            mark_segment = [interval = *result](const std::vector<double>& xs, const std::vector<double>& ys) {
                std::vector<double> out;
                for (size_t i = 0; i < xs.size(); ++i) {
                    if (xs[i] >= interval.left && xs[i] <= interval.right)
                        out.push_back(ys[i]);
                    else
                        out.push_back(std::numeric_limits<double>::quiet_NaN());
                }
                return out;
            };
        }

        if (result) {
            std::cout << "Interval of sign change: ["
                      << result->left << ", " << result->right << "]\n";
            std::cout << "f(left) = " << func(result->left)
                      << ", f(right) = " << func(result->right) << std::endl;
            std::cout << "Elapsed: " << (elapsed_us / 1000.0)
                      << " ms (" << elapsed_us << " us)\n";
        } else {
            std::cout << "No sign change found in interval ["
                      << a << ", " << b << "]\n";
        }

        if (plotter) {
            plotter->plot(
                x, y, "f(x) and interval", "x", "f(x)", false, mark_segment
            );
        }
    }

private:
    IRootFinder& finder;
    Plotter* plotter;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_ROOTFINDTASK_H