#ifndef NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEEXPLICIT_H
#define NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEEXPLICIT_H

#pragma once

#include "TaskWaveBase.h"

class TaskWaveExplicit : public TaskWaveBase {
public:
    explicit TaskWaveExplicit(Plotter* plotter = nullptr,
                              ExactFunc exact = nullptr)
        : TaskWaveBase(plotter, exact)
    {}

protected:
    const char* schemeName() const override {
        return "Explicit cross scheme";
    }

    void solveND(double tMax_nd,
                 double h_nd,
                 double tau_nd,
                 int N,
                 int M) override
    {
        int nx = N + 1;
        int nt = M + 1;

        double lambda = tau_nd / h_nd;
        if (lambda > 1.0) {
            std::cout << "[TaskWaveExplicit] Warning: CFL violated, lambda = "
                      << lambda << " > 1 (scheme may be unstable)\n";
        }

        Eigen::MatrixXd theta(nt, nx);
        theta.setZero();

        double xi0  = physParams.x0    / physParams.L;
        double d_xi = physParams.delta / physParams.L;

        double h_eff   = 1.0 / N;
        double tau_eff = tMax_nd / M;
        lambda = tau_eff / h_eff;

        for (int i = 0; i < nx; ++i) {
            double xi  = i * h_eff;
            double psi = 0.0;
            if (xi >= xi0 - d_xi && xi <= xi0 + d_xi)
                psi = 1.0;
            theta(1, i) = tau_eff * psi;
        }

        theta(0, 0)    = 0.0;
        theta(0, nx-1) = 0.0;
        theta(1, 0)    = 0.0;
        theta(1, nx-1) = 0.0;

        for (int s = 1; s < nt - 1; ++s) {
            for (int i = 1; i < nx - 1; ++i) {
                double yim1 = theta(s, i - 1);
                double yi   = theta(s, i);
                double yip1 = theta(s, i + 1);

                theta(s + 1, i) = 2.0 * yi - theta(s - 1, i)
                                  + lambda * lambda * (yip1 - 2.0 * yi + yim1);
            }
            theta(s + 1, 0)    = 0.0;
            theta(s + 1, nx-1) = 0.0;
        }

        double U0 = physParams.v0 * physParams.L / physParams.c;
        for (int s = 0; s < nt; ++s)
            for (int i = 0; i < nx; ++i)
                u(s, i) = U0 * theta(s, i);
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEEXPLICIT_H
