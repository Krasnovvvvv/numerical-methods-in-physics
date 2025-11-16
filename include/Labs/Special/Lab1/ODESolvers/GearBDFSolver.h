#ifndef NUMERICAL_METHODS_IN_PHYSICS_GEARBDFSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_GEARBDFSOLVER_H

#pragma once
#include "Base/IODESolver.h"
#include <vector>
#include <deque>
#include <string>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>

class GearBDFSolver : public IODESolver {
public:
    explicit GearBDFSolver(int order = 2, IODESolver* starter = nullptr, int max_iter = 20)
        : m_order(order), starter_solver(starter), max_iter(max_iter) {
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
        size_t N = y0.size();

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

            // Compute rhs of BDF equation (everything but y_{n+1})
            std::vector<double> rhs(N, 0.0);
            for(size_t j = 0; j < N; ++j) {
                rhs[j] = 0.0;
                for(size_t k = 1; k < a.size(); ++k)
                    rhs[j] += a[k] * y_hist[y_hist.size()-k][j];
                rhs[j] = -rhs[j];
            }

            // Newton's method for nonlinear system: F(y) = a0*y - h*b*f(t_{n+1}, y) - rhs = 0
            std::vector<double> y_new = state; // first guess (previous y)
            bool converged = false;
            for(int iter = 0; iter < max_iter; ++iter) {
                std::vector<double> fBDF = func(t_next, y_new);
                std::vector<double> F(N);
                for (size_t j = 0; j < N; ++j)
                    F[j] = a[0]*y_new[j] - h*b*fBDF[j] - rhs[j];

                // Check residual norm
                double norm = 0.0;
                for (double v : F) norm += std::abs(v);
                if(norm < tol) {
                    converged = true;
                    break;
                }

                // Finite-difference Jacobian
                std::vector<std::vector<double>> J(N, std::vector<double>(N, 0.0));
                double delta = 1e-8;
                for (size_t k = 0; k < N; ++k) {
                    std::vector<double> y_eps = y_new;
                    y_eps[k] += delta;
                    std::vector<double> f_eps = func(t_next, y_eps);
                    for (size_t j = 0; j < N; ++j)
                        J[j][k] = (a[0]*(y_eps[j]-y_new[j]) - h*b*(f_eps[j]-fBDF[j])) / delta;
                }

                // Solve J dx = -F
                std::vector<double> dx = solve_linear(J, F, -1.0);
                for(size_t j=0; j<N; ++j)
                    y_new[j] += dx[j];

                double dxnorm = 0.0;
                for(double x : dx) dxnorm += std::abs(x);
                if(dxnorm < tol) {
                    converged = true;
                    break;
                }
            }
            if(!converged)
                std::cerr << "(GearBDF) Warning: Newton's method did not converge at t=" << t_next << std::endl;

            // error estimation for logging/diagnostic
            double error = 0.0;
            for(size_t j = 0; j < N; ++j)
                error += std::abs(y_new[j] - state[j]);
            result.errorEstimates.push_back(error);

            state = y_new;
            time = t_next;
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

    // Universal n x n Gaussian elimination
    std::vector<double> solve_linear(std::vector<std::vector<double>> A, std::vector<double> b, double alpha = 1.0) {
        size_t n = b.size();
        for(size_t i=0; i<n; ++i) b[i] *= alpha;
        // Forward elimination
        for(size_t i=0; i<n; ++i) {
            size_t maxrow = i;
            for(size_t k=i+1; k<n; ++k)
                if(std::abs(A[k][i]) > std::abs(A[maxrow][i])) maxrow = k;
            std::swap(A[i], A[maxrow]);
            std::swap(b[i], b[maxrow]);
            for(size_t k=i+1; k<n; ++k) {
                double c = A[k][i]/A[i][i];
                for(size_t j=i; j<n; ++j)
                    A[k][j] -= c*A[i][j];
                b[k] -= c*b[i];
            }
        }
        // Back substitution
        std::vector<double> x(n);
        for(int i=int(n)-1; i>=0; --i) {
            x[i] = b[i]/A[i][i];
            for(int k=0; k<i; ++k)
                b[k] -= A[k][i]*x[i];
        }
        return x;
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_GEARBDFSOLVER_H
