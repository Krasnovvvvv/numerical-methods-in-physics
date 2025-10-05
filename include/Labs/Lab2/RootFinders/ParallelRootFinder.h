#ifndef NUMERICAL_METHODS_IN_PHYSICS_PARALLELROOTFINDER_H
#define NUMERICAL_METHODS_IN_PHYSICS_PARALLELROOTFINDER_H

#pragma once
#include "Base/IRootFinder.h"
#include <thread>
#include <vector>
#include <atomic>

class ParallelRootFinder final : public IRootFinder {
public:
    ParallelRootFinder(int numThreads = 4) : threads(numThreads) {}

    std::optional<RootInterval> findSignChange(
        std::function<double(double)> f,
        std::function<bool(double)> isInDomain,
        double a, double b, double step, double tol
    ) override {
        std::atomic<bool> found(false);
        std::optional<RootInterval> result = std::nullopt;
        double range = (b - a) / threads;
        std::vector<std::thread> pool;
        for (int t = 0; t < threads; ++t) {
            pool.emplace_back([&, t]() {
                double local_a = a + t * range;
                double local_b = (t == threads - 1) ? b : local_a + range;
                double x1 = local_a;
                double x2 = x1 + step;
                while (x2 <= local_b && !found.load()) {
                    if (!isInDomain(x1) || !isInDomain(x2)) {
                        x1 = x2;
                        x2 += step;
                        continue;
                    }
                    double y1 = f(x1), y2 = f(x2);
                    if ((y1 * y2 < -tol) || (std::abs(y1) < tol) || (std::abs(y2) < tol)) {
                        if (!found.exchange(true)) {
                            result = RootInterval{x1, x2};
                        }
                        return;
                    }
                    x1 = x2;
                    x2 += step;
                }
            });
        }
        for (auto &th : pool) th.join();
        return result;
    }
private:
    int threads;
};


#endif //NUMERICAL_METHODS_IN_PHYSICS_PARALLELROOTFINDER_H