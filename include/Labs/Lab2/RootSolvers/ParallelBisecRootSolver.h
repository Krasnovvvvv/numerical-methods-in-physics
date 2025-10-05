#ifndef NUMERICAL_METHODS_IN_PHYSICS_PARALLELROOTSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_PARALLELROOTSOLVER_H

#pragma once
#include "Base/IRootSolver.h"
#include "Helpers/refineRootBisection.h"
#include <thread>
#include <vector>
#include <atomic>
#include <optional>
#include <cmath>

class ParallelBisecRootSolver : public IRootSolver {
public:
    ParallelBisecRootSolver(int numThreads = 4) : threads(numThreads) {}

    std::optional<RootSolveResult> solve(
        std::function<double(double)> func,
        std::function<bool(double)> isInDomain,
        double tol,
        double x0,
        double x1,
        double step = 0.001,
        size_t max_iter = 100
    ) override {
        std::atomic<bool> found(false);
        std::optional<std::pair<double,double>> interval = std::nullopt;
        double range = (x1 - x0) / threads;
        std::vector<std::thread> pool;

        for (int t = 0; t < threads; ++t) {
            pool.emplace_back([&, t]() {
                double local_a = x0 + t * range;
                double local_b = (t == threads - 1) ? x1 : (local_a + range);
                double x1t = local_a;
                double x2t = x1t + step;
                while (x2t <= local_b && !found.load()) {
                    if (!isInDomain(x1t) || !isInDomain(x2t)) {
                        x1t = x2t;
                        x2t += step;
                        continue;
                    }
                    double y1 = func(x1t), y2 = func(x2t);
                    if ((y1 * y2 < -tol) || (std::abs(y1) < tol) || (std::abs(y2) < tol)) {
                        if (!found.exchange(true)) {
                            interval = std::make_pair(x1t, x2t);
                        }
                        return;
                    }
                    x1t = x2t;
                    x2t += step;
                }
            });
        }
        for (auto& th : pool) th.join();

        if (interval) {
            double root = refineRootBisection(func, interval->first, interval->second, tol, max_iter);
            return RootSolveResult{root, max_iter, std::abs(func(root))};
        }
        return std::nullopt;
    }
private:
    int threads;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_PARALLELROOTSOLVER_H
