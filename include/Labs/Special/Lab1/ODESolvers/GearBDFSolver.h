#ifndef NUMERICAL_METHODS_IN_PHYSICS_GEARBDFSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_GEARBDFSOLVER_H

#pragma once
#include "Base/IODESolver.h"
#include <vector>
#include <deque>
#include <string>
#include <cmath>
#include <stdexcept>

class GearBDFSolver : public IODESolver {
public:
    explicit GearBDFSolver(int order = 2, IODESolver* starter = nullptr, int max_iter = 50) :
    m_order(order), starter_solver(starter), max_iter(max_iter) {
        // Classic BDF coefficients
        if      (order == 1) { a = {1., -1.};           b = 1.; }
        else if (order == 2) { a = {1., -4./3., 1./3.};  b = 2./3.; }
        else if (order == 3) { a = {1., -18./11., 9./11., -2./11.}; b = 6./11.; }
        else if (order == 4) { a = {1., -48./25., 36./25., -16./25., 3./25.}; b = 12./25.; }
        else throw std::invalid_argument("GearBDF: Supported orders 1..4 only");
    }

    std::string name() const override { return "GearBDF (order " + std::to_string(m_order) + ")"; }

    ODEResult solve(
        const std::function<std::vector<double>(double, const std::vector<double>&)>& func,
        std::vector<double> y0, double t0, double tn, double h,
        double tol = 1e-8
    ) override {
        ODEResult result;
        double time = t0;
        std::vector<double> state = y0;
        std::deque<std::vector<double>> y_hist;

        // --- Initial fill of history using the starter solver ---
        int start_steps = m_order-1;
        if(start_steps > 0 && starter_solver) {
            ODEResult starter_result = starter_solver->solve(func, y0, t0, t0+h*start_steps, h, tol);
            for(size_t i = 0; i < starter_result.t.size(); ++i) {
                result.t.push_back(starter_result.t[i]);
                result.y.push_back(starter_result.y[i]);
                y_hist.push_back(starter_result.y[i]);
            }
            time = t0 + h*start_steps;
            state = starter_result.y.back();
        } else {
            for(int i = 0; i < start_steps && time <= tn + 1e-12; ++i) {
                result.t.push_back(time);
                result.y.push_back(state);
                y_hist.push_back(state);
                auto f = func(time, state);
                for(size_t j = 0; j < state.size(); ++j) state[j] += h*f[j];
                time += h;
            }
        }

        // --- Main BDF step ---
        while (time <= tn + 1e-12) {
            result.t.push_back(time);
            result.y.push_back(state);
            y_hist.push_back(state);
            if (y_hist.size() > m_order) y_hist.pop_front();

            double t_next = time + h;
            std::vector<double> y_next = state; // first guess

            // fixed-point iteration (f at t_{n+1}, y_{n+1})
            for(int iter = 0; iter < max_iter; ++iter) {
                auto f = func(t_next, y_next);
                std::vector<double> rhs(y0.size());
                for(size_t j = 0; j < rhs.size(); ++j) {
                    // sum: a[1]*y_n + a[2]*y_{n-1} + ...
                    rhs[j] = 0.0;
                    for(size_t k = 1; k < a.size(); ++k)
                        rhs[j] += a[k] * y_hist[y_hist.size()-k][j];
                    // add right part
                    rhs[j] = -rhs[j] + h*b*f[j];
                    rhs[j] /= a[0]; // divide by a0 (should be 1 for classic BDF)
                }
                // check for fixed-point convergence
                double diff = 0.0;
                for(size_t j = 0; j < y_next.size(); ++j)
                    diff += std::abs(rhs[j] - y_next[j]);
                y_next = rhs;
                if(diff < tol) break;
            }

            // error estimation
            double error = 0.0;
            for(size_t j = 0; j < state.size(); ++j)
                error += std::abs(y_next[j] - state[j]);
            result.errorEstimates.push_back(error);

            state = y_next;
            time += h;
        }
        result.steps = result.t.size();
        return result;
    }

private:
    int m_order;
    int max_iter;
    std::vector<double> a; // classical BDF coefficients
    double b; // right part coefficients
    IODESolver* starter_solver;
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_GEARBDFSOLVER_H
