#ifndef NUMERICAL_METHODS_IN_PHYSICS_TASKBESSEL_H
#define NUMERICAL_METHODS_IN_PHYSICS_TASKBESSEL_H
#pragma once

#include <Eigen/Dense>
#include "Base/ISolver.h"
#include "Helpers/Plotter.h"
#include "Helpers/Timer.h"
#include <cmath>
#include <iomanip>
#include <iostream>

class TaskBessel {
public:
    explicit TaskBessel(
        ISolver& solver,
        Plotter* plotter = nullptr
    )
    : solver(solver), plotter(plotter) {}

    void run(int order, double a, double b, double u_a,
        double u_b, int N
    ) {
        const size_t dim = N + 1;

        Eigen::VectorXd x(dim);
        Eigen::VectorXd u(dim);
        Eigen::MatrixXd A(dim, dim);
        Eigen::VectorXd B(dim);

        const double h = (b - a) / N;
        for (size_t i = 0; i <= N; ++i) {
            x[i] = a + i * h;
        }

        A.setZero(dim, dim);
        B.setZero(dim);

        A(0, 0) = 1.0;
        A(dim-1, dim-1) = 1.0;
        B[0] = u_a;
        B[dim - 1] = u_b;

        for (size_t i = 1; i < N; ++i) {
            double x_i = x[i];
            double x_ip1 = x[i + 1];
            double x_im1 = x[i - 1];
            double a_ip12 = std::pow((x_i + x_ip1) / 2.0, 2);
            double a_im12 = std::pow((x_i + x_im1) / 2.0, 2);
            double qi = std::pow(x_i, 2) - std::pow(order, 2);

            A(i, i - 1) = a_im12 / (h * h) + x_i / (2 * h);
            A(i, i) = - (a_im12 + a_ip12) / (h * h) + qi;
            A(i, i + 1) = a_ip12 / (h * h) - x_i / (2 * h);

            B[i] = 0.0;
        }

        Timer<std::chrono::microseconds> timer;
        auto result = solver.solve(A, B);
        auto elapsed_us = timer.elapsed();
        u = result.solution;

        std::vector<double> x_vec(dim), u_vec(dim);
        for (size_t i = 0; i <= N; ++i) {
            x_vec[i] = x[i];
            u_vec[i] = u[i];
        }

        std::cout << "--- The result of solving the Bessel " + std::to_string(order) + "order equation --- \n";
        std::cout << "Computation time: " << elapsed_us << " microseconds\n";
        std::cout << std::scientific;
        std::cout << "Residual: " << result.residual << "\n";

        if (plotter)
            plotter -> plot(x_vec, u_vec, "U(x) for " + std::to_string(order) +
                "D" + "Bessel equation", "x", "U");
    }
private:
    ISolver& solver;
    Plotter* plotter;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_TASKBESSEL_H