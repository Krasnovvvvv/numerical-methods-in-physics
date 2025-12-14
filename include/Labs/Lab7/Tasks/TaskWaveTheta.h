#ifndef NUMERICAL_METHODS_IN_PHYSICS_TASKWAVETHETA_H
#define NUMERICAL_METHODS_IN_PHYSICS_TASKWAVETHETA_H

#pragma once

#include "TaskWaveBase.h"
#include "Base/ISolver.h"

class TaskWaveTheta : public TaskWaveBase {
public:
    TaskWaveTheta(ISolver& solver,
                  Plotter* plotter = nullptr,
                  double sigma = 0.25,
                  ExactFunc exact = nullptr)
        : TaskWaveBase(plotter, exact)
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

        const double lambda2 = (tau_nd * tau_nd) / (h_nd * h_nd);

        Eigen::MatrixXd theta(nt, nx);
        theta.setZero();

        // --- начальные слои θ^0, θ^1 ---
        double xi0    = physParams.x0    / physParams.L;
        double d_xi   = physParams.delta / physParams.L;
        double h_eff  = 1.0 / N;
        double tau_eff = tMax_nd / M;

        // θ(ξ,0) = 0 уже стоит
        for (int i = 0; i < nx; ++i) {
            double xi  = i * h_eff;
            double psi = 0.0;
            if (xi >= xi0 - d_xi && xi <= xi0 + d_xi)
                psi = 1.0;                  // θ_τ = 1 на ударном отрезке
            theta(1, i) = tau_eff * psi;    // θ^1_i = τ * ψ_i
        }

        theta(0, 0)    = 0.0;
        theta(0, nx-1) = 0.0;
        theta(1, 0)    = 0.0;
        theta(1, nx-1) = 0.0;

        // --- матрица A для внутренних узлов (1..N-1) ---
        const int dim = nx - 2;
        Eigen::MatrixXd A(dim, dim);
        A.setZero();

        double a_off  = -sigma * lambda2;
        double a_diag =  2.0 + 2.0 * sigma * lambda2;

        for (int j = 0; j < dim; ++j) {
            if (j > 0)       A(j, j - 1) = a_off;
                              A(j, j)     = a_diag;
            if (j + 1 < dim) A(j, j + 1) = a_off;
        }

        // --- временной цикл ---
        for (int s = 1; s < nt - 1; ++s) {
            Eigen::VectorXd B(dim);
            B.setZero();

            for (int j = 0; j < dim; ++j) {
                int i = j + 1;  // внутренний узел

                double y_im1_s   = theta(s,     i - 1);
                double y_i_s     = theta(s,     i);
                double y_ip1_s   = theta(s,     i + 1);
                double y_im1_sm1 = theta(s - 1, i - 1);
                double y_i_sm1   = theta(s - 1, i);
                double y_ip1_sm1 = theta(s - 1, i + 1);

                double lap_s   =  y_ip1_s   - 2.0 * y_i_s   + y_im1_s;
                double lap_sm1 =  y_ip1_sm1 - 2.0 * y_i_sm1 + y_im1_sm1;

                // (1) классическая θ‑схема:
                // y^{s+1} - 2 y^s + y^{s-1} = λ² [ σ Λ y^{s+1}
                //                                + (1-2σ) Λ y^s
                //                                + σ Λ y^{s-1} ]
                //
                // после переноса σ Λ y^{s+1} влево:
                //
                //   (2 + 2σλ²) y_i^{s+1} + (-σλ²) y_{i-1}^{s+1} + (-σλ²) y_{i+1}^{s+1}
                //       = 2 y_i^s - y_i^{s-1}
                //         + λ² (1-2σ) Λ y_i^s + σλ² Λ y_i^{s-1}
                double rhs = 2.0 * y_i_s - y_i_sm1
                             + lambda2 * (1.0 - 2.0 * sigma) * lap_s
                             + sigma   * lambda2 * lap_sm1;

                B(j) = rhs;
            }

            auto res = solver.solve(A, B);

            for (int j = 0; j < dim; ++j) {
                int i = j + 1;
                theta(s + 1, i) = res.solution[j];
            }

            theta(s + 1, 0)    = 0.0;
            theta(s + 1, nx-1) = 0.0;
        }

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
