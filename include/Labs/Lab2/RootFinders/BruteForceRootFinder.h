#ifndef NUMERICAL_METHODS_IN_PHYSICS_BRUTEFORCEROOTFINDER_H
#define NUMERICAL_METHODS_IN_PHYSICS_BRUTEFORCEROOTFINDER_H

#pragma once
#include "Base/IRootFinder.h"

class BruteForceRootFinder final : public IRootFinder {
public:
    std::optional<RootInterval> findSignChange(
        std::function<double(double)> f,
        std::function<bool(double)> isInDomain,
        double a, double b, double step, double tol
    ) override {
        double x1 = a, x2 = a + step;
        while (x2 <= b) {
            if (!isInDomain(x1) || !isInDomain(x2)) {
                x1 = x2;
                x2 += step;
                continue;
            }
            double y1 = f(x1), y2 = f(x2);
            if ((y1 * y2 < -tol) || (std::abs(y1) < tol) || (std::abs(y2) < tol))
                return RootInterval{x1, x2};
            x1 = x2;
            x2 += step;
        }
        return std::nullopt;
    }
};



#endif //NUMERICAL_METHODS_IN_PHYSICS_BRUTEFORCEROOTFINDER_H