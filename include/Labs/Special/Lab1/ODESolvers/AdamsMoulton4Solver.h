#ifndef NUMERICAL_METHODS_IN_PHYSICS_ADAMSMOULTON4SOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_ADAMSMOULTON4SOLVER_H

#pragma once
#include "Base/IODESolver.h"
#include <deque>
#include <vector>

class AdamsMoulton4Solver : public IODESolver {
public:
    AdamsMoulton4Solver(IODESolver* starter_method = nullptr, int num_corrections = 1)
        : starter_solver(starter_method), corrections(num_corrections) {}

    std::string name() const override {
        return "Adams-Bashforth-Moulton 4 (" + std::to_string(corrections) + " corr.)";
    }

    ODEResult solve(
        const std::function<std::vector<double>(double, const std::vector<double>&)>& func,
        std::vector<double> y0, double t0, double tn, double h,
        double tol = 1e-8) override
    {
        ODEResult result;
        double time = t0;
        std::vector<double> state = y0;

        std::deque<std::vector<double>> f_hist;
        std::deque<std::vector<double>> y_hist;

        // --- Starting steps via delegated solver ---
        if (starter_solver) {
            ODEResult start = starter_solver->solve(func, y0, t0, t0 + 2.0 * h, h, tol);

            for (size_t i = 0; i < start.t.size(); ++i) {
                result.t.push_back(start.t[i]);
                result.y.push_back(start.y[i]);
                auto fi = func(start.t[i], start.y[i]);
                f_hist.push_front(fi);
                y_hist.push_front(start.y[i]);
            }
            time = t0 + 3.0 * h;
            state = start.y.back();
        } else {
            // Fallback: explicit Euler
            for (int i = 0; i < 3 && time <= tn; ++i) {
                result.t.emplace_back(time);
                result.y.emplace_back(state);
                auto f = func(time, state);
                f_hist.push_front(f);
                y_hist.push_front(state);

                for (size_t j = 0; j < state.size(); ++j)
                    state[j] += h * f[j];
                time += h;
            }
        }

        // --- The main cycle of the forecast-correction method ---
        while (time <= tn) {
            result.t.emplace_back(time);
            result.y.emplace_back(state);

            // The Adams-Bashfort predictor 4
            auto f0 = func(time, state);
            std::vector<double> predictor = state;
            for (size_t j = 0; j < state.size(); ++j) {
                double incr =
                    (55.0/24)*f0[j]
                    - (59.0/24)*f_hist[0][j]
                    + (37.0/24)*f_hist[1][j]
                    - (9.0/24)*f_hist[2][j];
                predictor[j] += h * incr;
            }

            // The Adams-Moulton Corrector 4
            std::vector<double> next_state = predictor;
            for (int iter = 0; iter < corrections; ++iter) {
                auto f_pred = func(time + h, next_state);
                std::vector<double> new_state = state;
                for (size_t j = 0; j < state.size(); ++j) {
                    double incr =
                        (9.0/24)*f_pred[j]
                        + (19.0/24)*f0[j]
                        - (5.0/24)*f_hist[0][j]
                        + (1.0/24)*f_hist[1][j];
                    new_state[j] += h * incr;
                }
                next_state.swap(new_state);
            }

            f_hist.push_front(f0);
            y_hist.push_front(state);
            while (f_hist.size() > 3) f_hist.pop_back();
            while (y_hist.size() > 3) y_hist.pop_back();

            state = next_state;
            time += h;
        }

        result.steps = result.t.size();
        return result;
    }

private:
    IODESolver* starter_solver;
    int corrections;
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_ADAMSMOULTON4SOLVER_H
