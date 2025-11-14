#ifndef NUMERICAL_METHODS_IN_PHYSICS_GEARSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_GEARSOLVER_H

#pragma once
#include "Base/IODESolver.h"
#include <deque>
#include <vector>
#include <stdexcept>

class GearSolver : public IODESolver {
public:
    GearSolver(int order = 2, IODESolver* starter_method = nullptr)
        : m_order(order), starter_solver(starter_method)
    {
        if (order < 1 || order > 4)
            throw std::invalid_argument("GearSolver: Only orders 1-4 supported.");
        gear_coeffs = {                // orders
            {1},                         // 1
            {3.0/2, -1.0/2},             // 2
            {23.0/12, -16.0/12, 5.0/12}, // 3
            {55.0/24, -59.0/24, 37.0/24, -9.0/24} // 4
        };
    }

    std::string name() const override { return "Gear"; }

    ODEResult solve(
        const std::function<std::vector<double>(double, const std::vector<double>&)>& func,
        std::vector<double> y0, double t0, double tn, double h,
        double tol = 1e-8) override {

        ODEResult result;
        int steps_req = m_order;
        double time = t0;
        std::vector<double> state = y0;

        std::deque<std::vector<double>> f_hist;
        std::deque<std::vector<double>> state_hist;

        // ---Starting steps via delegated solver ---
        if (starter_solver && steps_req > 1) {
            ODEResult start = starter_solver->solve(func, y0, t0, t0 + h * (steps_req - 2), h, tol);

            for (size_t i = 0; i < start.t.size(); ++i) {
                result.t.push_back(start.t[i]);
                result.y.push_back(start.y[i]);
                auto fi = func(start.t[i], start.y[i]);
                f_hist.push_front(fi);
                state_hist.push_front(start.y[i]);
            }
            time = t0 + h * (steps_req - 1);
            state = start.y.back();
        } else {
            // Fallback: explicit Euler
            for (int i = 0; i < steps_req-1 && time <= tn; ++i) {
                result.t.emplace_back(time);
                result.y.emplace_back(state);

                auto fi = func(time, state);
                f_hist.push_front(fi);
                state_hist.push_front(state);

                std::vector<double> prev_state = state;
                for (size_t j = 0; j < state.size(); ++j)
                    state[j] += h * fi[j];

                double error = 0;
                for (size_t j = 0; j < state.size(); ++j)
                    error += std::abs(state[j] - prev_state[j]);
                result.errorEstimates.emplace_back(error);

                time += h;
            }
        }

        // --- The main cycle of the Gear method ---
        while (time <= tn) {
            result.t.emplace_back(time);
            result.y.emplace_back(state);

            auto f0 = func(time, state);
            f_hist.push_front(f0);
            state_hist.push_front(state);
            while (f_hist.size() > m_order) f_hist.pop_back();
            while (state_hist.size() > m_order) state_hist.pop_back();

            std::vector<double> new_state = state_hist[0];
            for (size_t j = 0; j < new_state.size(); ++j) {
                double incr = 0;
                for (int k = 0; k < m_order; ++k)
                    incr += gear_coeffs[m_order-1][k] * f_hist[k][j];
                new_state[j] += h * incr;
            }

            double error = 0;
            for (size_t j = 0; j < new_state.size(); ++j)
                error += std::abs(new_state[j] - state[j]);
            result.errorEstimates.emplace_back(error);

            state = new_state;
            time += h;
        }

        result.steps = result.t.size();
        return result;
    }

private:
    int m_order;
    std::vector<std::vector<double>> gear_coeffs;
    IODESolver* starter_solver;
};
#endif // NUMERICAL_METHODS_IN_PHYSICS_GEARSOLVER_H
