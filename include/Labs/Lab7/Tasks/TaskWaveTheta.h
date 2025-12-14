#ifndef NUMERICAL_METHODS_IN_PHYSICS_TASKWAVETHETA_H
#define NUMERICAL_METHODS_IN_PHYSICS_TASKWAVETHETA_H

#pragma once

#include "TaskWaveBase.h"
#include "Base/ISolver.h"

class TaskWaveTheta : public TaskWaveBase {
public:
    // 0 <= sigma <= 1/2; при sigma >= 1/4 схема безусловно устойчива
    TaskWaveTheta(ISolver& solver,
                  Plotter* plotter = nullptr,
                  double sigma = 0.25)
        : TaskWaveBase(plotter)
        , solver(solver)
        , sigma(sigma)
    {}

protected:
    const char* schemeName() const override {
        return "Theta-scheme for wave (weights)";
    }

    void solveND(double tMax_nd,
                 double h_nd,
                 double tau_nd,
                 int N,
                 int M) override
    {
        const int nx = N + 1;
        const int nt = M + 1;

        const double lambda2 = (tau_nd * tau_nd) / (h_nd * h_nd); // (c τ / h)^2, c=1 в ND
        const double a = sigma * lambda2;

        Eigen::MatrixXd theta(nt, nx);
        theta.setZero();

        double xi0  = physParams.x0    / physParams.L;
        double d_xi = physParams.delta / physParams.L;
        double h_eff   = 1.0 / N;
        double tau_eff = tMax_nd / M;

        for (int i = 0; i < nx; ++i) {
            double xi  = i * h_eff;
            double psi = 0.0;
            if (xi >= xi0 - d_xi && xi <= xi0 + d_xi)
                psi = 1.0;
            theta(1, i) = tau_eff * psi;    // θ^1_i
        }
        theta(0, 0)    = 0.0;
        theta(0, nx-1) = 0.0;
        theta(1, 0)    = 0.0;
        theta(1, nx-1) = 0.0;

        const int dim = nx - 2;
        Eigen::MatrixXd A(dim, dim);
        A.setZero();

        for (int j = 0; j < dim; ++j) {
            if (j > 0)         A(j, j-1) = -a;
                               A(j, j)   = 2.0 + 2.0 * a;
            if (j + 1 < dim)   A(j, j+1) = -a;
        }

        // временной цикл: на каждом шаге решаем A * y^{s+1}_int = RHS
        for (int s = 1; s < nt - 1; ++s) {
            Eigen::VectorXd B(dim);
            B.setZero();

            for (int j = 0; j < dim; ++j) {
                int i = j + 1; // внутренний индекс 1..N-1

                double y_im1_s = theta(s, i - 1);
                double y_i_s   = theta(s, i);
                double y_ip1_s = theta(s, i + 1);
                double y_i_sm1 = theta(s - 1, i);

                // Λ y^s_i = (y_{i+1}^s - 2 y_i^s + y_{i-1}^s) / h^2,
                double lap_num = y_ip1_s - 2.0 * y_i_s + y_im1_s;

                double rhs = 2.0 * y_i_s - y_i_sm1
                             + lambda2 * (1.0 - 2.0 * sigma) * lap_num;

                B(j) = rhs;
            }

            auto res = solver.solve(A, B);

            // обновляем внутренние узлы слоя s+1
            for (int j = 0; j < dim; ++j) {
                int i = j + 1;
                theta(s + 1, i) = res.solution[j];
            }
            // граничные условия Дирихле
            theta(s + 1, 0)    = 0.0;
            theta(s + 1, nx-1) = 0.0;
        }

        // перевод θ -> u
        double U0 = physParams.v0 * physParams.L / physParams.c;
        for (int s = 0; s < nt; ++s)
            for (int i = 0; i < nx; ++i)
                u(s, i) = U0 * theta(s, i);
    }

private:
    ISolver& solver;
    double sigma;
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_TASKWAVETHETA_H
