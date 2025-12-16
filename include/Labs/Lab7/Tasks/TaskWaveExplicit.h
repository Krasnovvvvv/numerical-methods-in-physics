#ifndef NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEEXPLICIT_H
#define NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEEXPLICIT_H

#pragma once

#include "TaskWaveBase.h"
#include <cmath>
#include <iostream>

class TaskWaveExplicit : public TaskWaveBase {
public:
    explicit TaskWaveExplicit(Plotter* plotter = nullptr,
                             ExactFunc exact = nullptr)
        : TaskWaveBase(plotter, exact)
    {}

protected:
    const char* schemeName() const override {
        return "Explicit cross scheme (sigma=0)";
    }

    void solveND(double /*tMax*/,
                 double h_x,
                 double tau_t,
                 int N,
                 int M) override
    {
        const int nx = N + 1;
        const int nt = M + 1;

        // lambda = c * tau / h
        const double lambda = physParams.c * tau_t / h_x;
        const double lambda2 = lambda * lambda;

        if (lambda > 1.0) {
            std::cout << "[TaskWaveExplicit] WARNING: CFL violated, lambda = "
                      << lambda << " > 1\n";
        }

        Eigen::MatrixXd y(nt, nx);
        y.setZero();

        // ---- Начальные условия ----
        // u(x,0) = 0
        // u_t(x,0) = v0 * chi_{[x0-d, x0+d]}(x)
        //
        // На сетке первый слой по времени:
        // y^1_i = tau * u_t(x_i, 0) = tau * v0 * chi_i

        for (int i = 0; i < nx; ++i) {
            double x_i = i * h_x;
            double chi = 0.0;
            if (x_i >= physParams.x0 - physParams.delta &&
                x_i <= physParams.x0 + physParams.delta)
                chi = 1.0;

            y(1, i) = tau_t * physParams.v0 * chi;
        }

        // Граничные условия Дирихле: y = 0 на концах
        y(0, 0)      = 0.0;
        y(0, nx - 1) = 0.0;
        y(1, 0)      = 0.0;
        y(1, nx - 1) = 0.0;

        // ---- Явная схема "крест" ----
        // y^{s+1}_i = lambda^2 * (y^s_{i+1} - 2*y^s_i + y^s_{i-1})
        //             + 2*y^s_i - y^{s-1}_i

        for (int s = 1; s < nt - 1; ++s) {
            for (int i = 1; i < nx - 1; ++i) {
                double y_sp = y(s,     i);
                double y_sm = y(s - 1, i);
                double y_ip = y(s, i + 1);
                double y_im = y(s, i - 1);

                y(s + 1, i) = lambda2 * (y_ip - 2.0 * y_sp + y_im)
                            + 2.0 * y_sp - y_sm;
            }

            // Граничные условия
            y(s + 1, 0)      = 0.0;
            y(s + 1, nx - 1) = 0.0;
        }

        u = y;
    }
};

#endif
