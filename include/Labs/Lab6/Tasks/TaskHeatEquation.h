#ifndef NUMERICAL_METHODS_IN_PHYSICS_TASKHEATEQ_H
#define NUMERICAL_METHODS_IN_PHYSICS_TASKHEATEQ_H

#pragma once

#include "Base/ISolver.h"
#include "Helpers/Plotter.h"
#include "Helpers/Timer.h"

#include <Eigen/Dense>
#include <vector>
#include <functional>
#include <cmath>
#include <iostream>

class TaskHeatEquation {
public:
    using FuncXT = std::function<double(double,double)>;
    using FuncX  = std::function<double(double)>;
    using FuncT  = std::function<double(double)>;

    explicit TaskHeatEquation(
        ISolver& solver,
        Plotter* plotter = nullptr,
        double sigma = 0.0           // 0 – явная (FTCS), 0.5 – Кранк–Николсон, 1 – неявная
    )
        : solver(solver), plotter(plotter), sigma(sigma) {}

    void run(
        double kappa,
        double l, double tMax,
        int N, int M,
        FuncXT f,
        FuncX u0,
        FuncT nu1, FuncT nu2,
        FuncXT exactSolution = nullptr
    )
    {
        const size_t nx = static_cast<size_t>(N) + 1;
        const size_t nt = static_cast<size_t>(M) + 1;

        const double h   = l / N;
        const double tau = tMax / M;

        // для уравнения u_t = -kappa * u_xx + f
        const double mu  = kappa * tau / (h * h);   // κ τ / h²

        Eigen::VectorXd x(nx), t(nt);
        for (size_t i = 0; i < nx; ++i) x[i] = i * h;
        for (size_t s = 0; s < nt; ++s) t[s] = s * tau;

        Eigen::MatrixXd u(nt, nx);

        // начальное условие
        for (size_t i = 0; i < nx; ++i)
            u(0, i) = u0(x[i]);

        Timer timer;

        for (size_t s = 0; s < nt - 1; ++s) {
            double ts   = t[s];
            double tsp1 = t[s+1];

            // граничные условия (Дирихле) на шаге n и n+1
            u(s,   0)    = nu1(ts);
            u(s,   nx-1) = nu2(ts);
            u(s+1, 0)    = nu1(tsp1);
            u(s+1, nx-1) = nu2(tsp1);

            if (sigma == 0.0) {
                // явная схема для u_t = -kappa u_xx + f:
                // u_i^{n+1} = u_i^n - μ (u_{i+1}^n - 2u_i^n + u_{i-1}^n) + τ f_i^n
                // нужно выбирать μ <= 0.5
                for (size_t i = 1; i < nx-1; ++i) {
                    double ui   = u(s,i);
                    double uim1 = u(s,i-1);
                    double uip1 = u(s,i+1);
                    double lap  = uip1 - 2.0*ui + uim1;
                    u(s+1,i) = ui + mu * lap + tau * f(x[i], ts);
                }
            } else {
                // θ‑схема (sigma=θ) для u_t = -kappa u_xx + f:
                // (1 + 2θμ)u_i^{n+1} - θμ u_{i-1}^{n+1} - θμ u_{i+1}^{n+1}
                // = (1 - 2(1-θ)μ)u_i^{n} + (1-θ)μ(u_{i-1}^n + u_{i+1}^n) + τ f_i^n
                const size_t dim = nx - 2;       // внутренние узлы
                Eigen::MatrixXd A(dim, dim);
                Eigen::VectorXd B(dim);
                A.setZero(); B.setZero();

                for (size_t j = 0; j < dim; ++j) {
                    size_t i = j + 1;           // индекс по x (1..nx-2)

                    // матрица для шага n+1
                    if (j > 0)        A(j, j-1) = -sigma * mu;
                    A(j, j) = 1.0 + 2.0 * sigma * mu;
                    if (j + 1 < dim)  A(j, j+1) = -sigma * mu;

                    // правая часть
                    double uim1 = u(s,i-1);
                    double ui   = u(s,i);
                    double uip1 = u(s,i+1);
                    double lap_n = uip1 - 2.0*ui + uim1;

                    // знак минус перед лапласианом
                    double rhs = ui + (1.0 - sigma) * mu * lap_n + tau * f(x[i], ts);

                    // учёт граничных значений на шаге n+1
                    if (i == 1) {
                        rhs += sigma * mu * u(s+1, 0);         // сосед слева
                    }
                    if (i == nx-2) {
                        rhs += sigma * mu * u(s+1, nx-1);      // сосед справа
                    }

                    B(j) = rhs;
                }

                auto result = solver.solve(A, B);
                for (size_t j = 0; j < dim; ++j)
                    u(s+1, j+1) = result.solution[j];
            }
        }

        auto elapsed_us = timer.elapsed();

        std::cout << "--- Heat equation solution ---\n";
        std::cout << "sigma = " << sigma
                  << ", mu = " << mu
                  << ", computation time: " << elapsed_us << " microseconds\n";

        if (plotter) {
            // выбираем более «удачные» моменты времени (доли от tMax)
            std::vector<double> rel_times = {0.1};

            // общий X для всех кривых
            std::vector<double> x_vec(nx);
            for (size_t i = 0; i < nx; ++i)
                x_vec[i] = x[i];

            std::vector<std::vector<double>> xs;
            std::vector<std::vector<double>> ys;
            std::vector<std::string>          labels;

            for (double r : rel_times) {
                double target_t = r * tMax;
                size_t s = static_cast<size_t>(std::round(target_t / tau));
                if (s >= nt) s = nt - 1;
                double ts = t[s];

                std::vector<double> u_num(nx);
                for (size_t i = 0; i < nx; ++i)
                    u_num[i] = u(s, i);

                xs.push_back(x_vec);
                ys.push_back(u_num);
                labels.emplace_back("u_{num}(x,t=" + std::to_string(ts) + ")");

                if (exactSolution) {
                    std::vector<double> u_ex(nx);
                    for (size_t i = 0; i < nx; ++i)
                        u_ex[i] = exactSolution(x[i], ts);

                    xs.push_back(x_vec);
                    ys.push_back(u_ex);
                    labels.emplace_back("u_{exact}(x,t=" + std::to_string(ts) + ")");
                }
            }

            plotter->plot(xs, ys, labels, "x", "u");
        }
    }

private:
    ISolver& solver;
    Plotter* plotter;
    double sigma;
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_TASKHEATEQ_H
