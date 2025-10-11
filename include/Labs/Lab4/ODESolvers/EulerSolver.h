#ifndef NUMERICAL_METHODS_IN_PHYSICS_EULERSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_EULERSOLVER_H

#pragma once
#include "Base/IODESolver.h"

class EulerSolver : public IODESolver {
public:
    std::string name() const override { return "Euler"; }

    ODEResult solve(
    const std::function<std::vector<double>(double, const std::vector<double>&)>& func,
    std::vector<double> y0, double t0, double tn, double h,
    double tol = 1e-8) override {
        ODEResult result;
        double time = t0;
        std::vector<double> state = y0;
        std::vector<double> prev_state = y0;

        while (time <= tn) {
            result.t.emplace_back(time);
            result.y.emplace_back(state);

            auto k1 = func(time, state);
            prev_state = state;
            for (size_t i = 0; i < state.size(); ++i)
                state[i] += h * k1[i];

            double error = 0;
            for (size_t i = 0; i < state.size(); ++i)
                error += std::abs(state[i] - prev_state[i]);
            result.errorEstimates.emplace_back(error);

            if (error < tol)
                break;

            time += h;
        }
        result.steps = result.t.size();
        return result;
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_EULERSOLVER_H