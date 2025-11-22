#ifndef NUMERICAL_METHODS_IN_PHYSICS_HEATRODTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_HEATRODTASK_H
#pragma once

#include <Eigen/Dense>
#include "Base/ISolver.h"
#include "Helpers/Plotter.h"
#include "Helpers/Timer.h"
#include <cmath>
#include <iomanip>
#include <iostream>

class HeatRodTask {
public:
    explicit HeatRodTask(
        ISolver& solver,
        Plotter* plotter = nullptr
    )
    : solver(solver), plotter(plotter) {}

    void run(
        const std::vector<double> &hP_over_kA,
        double L,
        const std::vector<int>& Ns,
        double Tb,
        double T_inf
    ) {
        using std::vector;
        using std::string;
        vector<vector<double>> xs, ys;
        vector<string> labels;

        for (double coeff : hP_over_kA) {
            double Bi = coeff * L * L;

            for (int N : Ns) {
                const size_t dim = N + 1;
                Eigen::VectorXd xi(dim);
                Eigen::VectorXd theta(dim);
                Eigen::MatrixXd A(dim, dim);
                Eigen::VectorXd B(dim);

                const double h = 1.0 / N;
                for (size_t i = 0; i <= N; ++i)
                    xi[i] = i * h;

                A.setZero(dim, dim);
                B.setZero(dim);

                A(0, 0) = 1.0;
                B[0] = 1.0;

                for (size_t i = 1; i < N; ++i) {
                    A(i, i - 1) = 1.0 / (h * h);
                    A(i, i) = -2.0 / (h * h) - Bi;
                    A(i, i + 1) = 1.0 / (h * h);
                    B[i] = 0.0;
                }

                A(dim - 1, dim - 1) = 1.0;
                A(dim - 1, dim - 2) = -1.0;
                B[dim - 1] = 0.0;

                Timer<std::chrono::microseconds> timer;
                auto result = solver.solve(A, B);
                auto elapsed_us = timer.elapsed();
                theta = result.solution;

                vector<double> x_vec(dim), T_vec(dim);
                for (size_t i = 0; i < dim; ++i) {
                    x_vec[i] = xi[i] * L;
                    T_vec[i] = T_inf + (Tb - T_inf) * theta[i];
                }

                xs.push_back(std::move(x_vec));
                ys.push_back(std::move(T_vec));
                labels.push_back("N = " + std::to_string(N) + ", Bi = " + std::to_string(Bi));

                std::cout << "N = " << N
                          << " | Time: " << elapsed_us << " microseconds"
                          << " | Residual: " << result.residual << "\n";
            }
        }

        if (plotter) {
            plotter->plot(xs, ys, labels, "x, m", "T(x), K");
        }
    }

private:
    ISolver& solver;
    Plotter* plotter;
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_HEATRODTASK_H