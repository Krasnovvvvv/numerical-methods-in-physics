#ifndef NUMERICAL_METHODS_IN_PHYSICS_BRUTEFORCEROOTFINDER_H
#define NUMERICAL_METHODS_IN_PHYSICS_BRUTEFORCEROOTFINDER_H

#pragma once
#include "Base/IRootSolver.h"
#include "Helpers/refineRootBisection.h"

class BruteForceBisecRootSolver : public IRootSolver {
public:
    std::optional<RootSolveResult> solve(
        std::function<double(double)> func,
        std::function<bool(double)> isInDomain,
        double tol,
        double x0,
        double x1 = 1,
        double step = 0.001,
        size_t max_iter = 100
    ) override {
        size_t iter = 0;
        x1 = x0 + step;
        while (x0 <= x1) {
            if (!isInDomain(x0) || !isInDomain(x1)) {
                x0 = x1;
                x1 += step;
                continue;
            }
            double y1 = func(x0), y2 = func(x1);
            if ((y1 * y2 < -tol) || (std::abs(y1) < tol) || (std::abs(y2) < tol)) {
                auto root = refineRootBisection(func,x0,x1,tol,max_iter);
                return RootSolveResult{root, iter,std::abs(func(root))};
            }
            x0 = x1;
            x1 += step;
            ++iter;
        }
        return std::nullopt;
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_BRUTEFORCEROOTFINDER_H