#ifndef NUMERICAL_METHODS_IN_PHYSICS_IMPROVEDEULERSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_IMPROVEDEULERSOLVER_H

#pragma once
#include "Base/IODESolver.h"

class ImprovedEulerSolver : public IODESolver {
public:
    std::string name() const override { return "Improved Euler"; }

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
            std::vector<double> predictor(state.size());
            for (size_t i = 0; i < state.size(); ++i)
                predictor[i] = state[i] + h * k1[i];
            auto k2 = func(time + h, predictor);

            prev_state = state;
            for (size_t i = 0; i < state.size(); ++i)
                state[i] += h * 0.5 * (k1[i] + k2[i]);

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

#endif //NUMERICAL_METHODS_IN_PHYSICS_IMPROVEDEULERSOLVER_H