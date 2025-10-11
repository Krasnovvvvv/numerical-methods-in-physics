#ifndef NUMERICAL_METHODS_IN_PHYSICS_RUNGEKUTTA4SOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_RUNGEKUTTA4SOLVER_H

#pragma once
#include "Base/IODESolver.h"

class RungeKutta4Solver : public IODESolver {
public:
    std::string name() const override { return "Runge-Kutta 4"; }

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

            std::vector<double> yk2(state.size());
            for (size_t i = 0; i < state.size(); ++i)
                yk2[i] = state[i] + h * 0.5 * k1[i];
            auto k2 = func(time + 0.5 * h, yk2);

            std::vector<double> yk3(state.size());
            for (size_t i = 0; i < state.size(); ++i)
                yk3[i] = state[i] + h * 0.5 * k2[i];
            auto k3 = func(time + 0.5 * h, yk3);

            std::vector<double> yk4(state.size());
            for (size_t i = 0; i < state.size(); ++i)
                yk4[i] = state[i] + h * k3[i];
            auto k4 = func(time + h, yk4);

            prev_state = state;
            for (size_t i = 0; i < state.size(); ++i)
                state[i] += h / 6.0 * (k1[i] +
                    2 * k2[i] + 2 * k3[i] + k4[i]);

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

#endif //NUMERICAL_METHODS_IN_PHYSICS_RUNGEKUTTA4SOLVER_H