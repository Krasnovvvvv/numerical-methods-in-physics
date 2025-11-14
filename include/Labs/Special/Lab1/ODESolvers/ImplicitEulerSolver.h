#ifndef NUMERICAL_METHODS_IN_PHYSICS_IMPLICITEULERSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_IMPLICITEULERSOLVER_H

#pragma once
#include "Base/IODESolver.h"
#include <vector>
#include <string>
#include <cmath>

class ImplicitEulerSolver : public IODESolver {
public:
    std::string name() const override { return "Implicit Euler"; }

    ODEResult solve(
        const std::function<std::vector<double>(double, const std::vector<double>&)>& func,
        std::vector<double> y0, double t0, double tn, double h,
        double tol = 1e-8, int max_iter = 50
    ) override {
        ODEResult result;
        double time = t0;
        std::vector<double> state = y0;

        while (time <= tn) {
            result.t.emplace_back(time);
            result.y.emplace_back(state);

            double next_time = time + h;
            // First approximation: explicit Euler step
            std::vector<double> y_next = state;
            auto k0 = func(time, state);
            for (size_t i = 0; i < y_next.size(); ++i)
                y_next[i] += h * k0[i];

            // the simplest iterations for searching y_{n+1} = y_n + h * f(t_{n+1}, y_{n+1})
            for (int iter = 0; iter < max_iter; ++iter) {
                auto k1 = func(next_time, y_next);
                std::vector<double> y_new = state;
                for (size_t i = 0; i < y_new.size(); ++i)
                    y_new[i] += h * k1[i];
                double diff = 0;
                for (size_t i = 0; i < y_new.size(); ++i)
                    diff += std::abs(y_new[i] - y_next[i]);
                y_next = y_new;
                if (diff < tol) break;
            }

            double error = 0.0;
            for (size_t i = 0; i < y_next.size(); ++i)
                error += std::abs(y_next[i] - state[i]);
            result.errorEstimates.emplace_back(error);

            state = y_next;
            time = next_time;
        }

        result.steps = result.t.size();
        return result;
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_IMPLICITEULERSOLVER_H