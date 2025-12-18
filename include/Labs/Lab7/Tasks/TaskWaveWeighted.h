#ifndef NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEWEIGHTED_H
#define NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEWEIGHTED_H

#pragma once

#include "TaskWaveBase.h"
#include <Eigen/Dense>
#include <iostream>

class TaskWaveWeighted : public TaskWaveBase {
public:
    explicit TaskWaveWeighted(Plotter* plotter = nullptr,
                              ExactFunc exact = nullptr,
                              PlotMode mode = PLOT_SOLUTION,
                              double sigma = 0.25)   // вес σ
        : TaskWaveBase(plotter, exact, mode)
        , sigma_(sigma)
    {}

protected:
    const char* schemeName() const override {
        return "Weighted scheme (sigma)";
    }

    void solveND(double /*tMax*/,
                 double h_x,
                 double tau_t,
                 int N,
                 int M) override
    {
        const int nx = N + 1;
        const int nt = M + 1;

        const double c  = physParams.c;
        const double lambda2 = (c * c) * (tau_t * tau_t) / (h_x * h_x);

        Eigen::MatrixXd y(nt, nx);
        y.setZero();

        // ---- начальные условия, как в явной схеме ----
        for (int i = 0; i < nx; ++i) {
            double x_i = i * h_x;
            double chi = 0.0;
            if (x_i >= physParams.x0 - physParams.delta &&
                x_i <= physParams.x0 + physParams.delta)
                chi = 1.0;
            y(1, i) = tau_t * physParams.v0 * chi;
        }

        // ГУ Дирихле
        y(0, 0)      = 0.0; y(0, nx - 1) = 0.0;
        y(1, 0)      = 0.0; y(1, nx - 1) = 0.0;

        // Коэффициенты для трёхдиагональной матрицы
        const double a = -sigma_ * lambda2;            // поддиагональ и наддиагональ
        const double b = 1.0 + 2.0 * sigma_ * lambda2; // диагональ

        Eigen::VectorXd rhs(nx - 2);
        Eigen::VectorXd y_new(nx - 2);

        // ---- шаги по времени ----
        for (int s = 1; s < nt - 1; ++s) {
            // правая часть для внутренних узлов i = 1..nx-2
            for (int i = 1; i < nx - 1; ++i) {
                double y_i_s   = y(s, i);
                double y_ip_s  = y(s, i + 1);
                double y_im_s  = y(s, i - 1);
                double y_i_sm1 = y(s - 1, i);
                double y_ip_sm1 = y(s - 1, i + 1);
                double y_im_sm1 = y(s - 1, i - 1);

                double term_center = 2.0 * y_i_s - y_i_sm1;
                double term_space  =
                    lambda2 * ((1.0 - 2.0 * sigma_) * (y_ip_s - 2.0 * y_i_s + y_im_s) +
                               sigma_ * (y_ip_sm1 - 2.0 * y_i_sm1 + y_im_sm1));

                rhs[i - 1] = term_center + term_space;
            }

            // учёт граничных узлов в правой части
            // (a*y_{0}^{s+1} и a*y_{nx-1}^{s+1} занулены, т.к. ГУ Дирихле: 0)

            // решение трёхдиагональной СЛАУ (метод прогонки)
            triDiagSolve(a, b, a, rhs, y_new);

            // записываем новый слой
            y(s + 1, 0)      = 0.0;
            y(s + 1, nx - 1) = 0.0;
            for (int i = 1; i < nx - 1; ++i)
                y(s + 1, i) = y_new[i - 1];
        }

        u = y;
    }

private:
    double sigma_;

    // решение трёхдиагональной системы Ax = d,
    // A: diag=b, subdiag=a, supdiag=c (здесь a=c)
    static void triDiagSolve(double a, double b, double c,
                             const Eigen::VectorXd& d,
                             Eigen::VectorXd& x)
    {
        int n = static_cast<int>(d.size());
        x.resize(n);

        Eigen::VectorXd alpha(n), beta(n);

        // прямой ход
        alpha[0] = -c / b;
        beta[0]  = d[0] / b;
        for (int i = 1; i < n; ++i) {
            double denom = b + a * alpha[i - 1];
            alpha[i] = -c / denom;
            beta[i]  = (d[i] - a * beta[i - 1]) / denom;
        }

        // обратный ход
        x[n - 1] = beta[n - 1];
        for (int i = n - 2; i >= 0; --i)
            x[i] = alpha[i] * x[i + 1] + beta[i];
    }
};

#endif
