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

    enum class PlotMode {
        Distribution = 1, // график распределения
        ErrorStudy   = 2  // исследование ошибки
    };

    explicit TaskHeatEquation(
        ISolver& solver,
        Plotter* plotter = nullptr,
        double sigma = 0.0,
        PlotMode mode = PlotMode::Distribution
    )
        : solver(solver)
        , plotter(plotter)
        , sigma(sigma)
        , mode(mode)
    {}

    void run(
        double kappa,
        double l, double tMax,
        double h, double tau,
        FuncXT f,
        FuncX u0,
        FuncT nu1, FuncT nu2,
        FuncXT exactSolution = nullptr
    ) {
        const int N = static_cast<int>(std::round(l / h));
        const int M = static_cast<int>(std::round(tMax / tau));
        runGrid(kappa, l, tMax, N, M, f, u0, nu1, nu2, exactSolution);
    }

private:
    void runGrid(
        double kappa,
        double l, double tMax,
        int N, int M,
        FuncXT f,
        FuncX u0,
        FuncT nu1, FuncT nu2,
        FuncXT exactSolution
    ) {
        const size_t nx  = static_cast<size_t>(N) + 1;
        const size_t nt  = static_cast<size_t>(M) + 1;
        const double h   = l / N;
        const double tau = tMax / M;
        const double mu  = kappa * tau / (h * h);

        Eigen::VectorXd x(nx), t(nt);
        for (size_t i = 0; i < nx; ++i)  x[i] = i * h;
        for (size_t s = 0; s < nt; ++s)  t[s] = s * tau;

        Eigen::MatrixXd u(nt, nx);
        for (size_t i = 0; i < nx; ++i)
            u(0, i) = u0(x[i]);

        Timer timer;

        for (size_t s = 0; s < nt - 1; ++s) {
            double ts   = t[s];
            double tsp1 = t[s + 1];

            u(s, 0)      = nu1(ts);
            u(s, nx - 1) = nu2(ts);
            u(s + 1, 0)      = nu1(tsp1);
            u(s + 1, nx - 1) = nu2(tsp1);

            if (sigma == 0.0) {
                for (size_t i = 1; i < nx - 1; ++i) {
                    double ui   = u(s, i);
                    double uim1 = u(s, i - 1);
                    double uip1 = u(s, i + 1);
                    double lap  = uip1 - 2.0 * ui + uim1;
                    u(s + 1, i) = ui + mu * lap + tau * f(x[i], ts);
                }
            } else {
                const size_t dim = nx - 2;
                Eigen::MatrixXd A(dim, dim);
                Eigen::VectorXd B(dim);
                A.setZero(); B.setZero();

                for (size_t j = 0; j < dim; ++j) {
                    size_t i = j + 1;

                    if (j > 0)       A(j, j - 1) = -sigma * mu;
                    A(j, j) = 1.0 + 2.0 * sigma * mu;
                    if (j + 1 < dim) A(j, j + 1) = -sigma * mu;

                    double uim1 = u(s, i - 1);
                    double ui   = u(s, i);
                    double uip1 = u(s, i + 1);
                    double lap_n = uip1 - 2.0 * ui + uim1;

                    double rhs = ui + (1.0 - sigma) * mu * lap_n
                                 + tau * f(x[i], ts);

                    if (i == 1)
                        rhs += sigma * mu * u(s + 1, 0);
                    if (i == nx - 2)
                        rhs += sigma * mu * u(s + 1, nx - 1);

                    B(j) = rhs;
                }

                auto result = solver.solve(A, B);
                for (size_t j = 0; j < dim; ++j)
                    u(s + 1, j + 1) = result.solution[j];
            }
        }

        auto elapsed_us = timer.elapsed();
        std::cout << "--- Heat equation solution ---\n";
        std::cout << "sigma = " << sigma
                  << ", mu = " << kappa * (t[1] - t[0]) / (h * h)
                  << ", computation time: " << elapsed_us << " microseconds\n";

        if (!plotter) return;

        if (mode == PlotMode::Distribution) {
            plotDistribution(u, x, t, tMax, exactSolution);
        } else if (mode == PlotMode::ErrorStudy && exactSolution) {
            plotDistribution(u, x, t, tMax, exactSolution);
        }
    }

    void plotDistribution(
        const Eigen::MatrixXd& u,
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& t,
        double tMax,
        FuncXT exactSolution
    ) {
        const size_t nx = static_cast<size_t>(x.size());
        const size_t nt = static_cast<size_t>(t.size());
        double tau = t[1] - t[0];

        std::vector<double> rel_times = {0.5, 1.0};

        std::vector<double> x_vec(nx);
        for (size_t i = 0; i < nx; ++i)
            x_vec[i] = x[i];

        std::vector<std::vector<double>> xs;
        std::vector<std::vector<double>> ys;
        std::vector<std::string> labels;

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

public:
    // эксперимент: ошибка по шагу h и по σ
    void runErrorStudy(
        double kappa,
        double l, double tMax,
        double h0, double tau0,
        FuncXT f,
        FuncX u0,
        FuncT nu1, FuncT nu2,
        FuncXT exactSolution
    ) {
        if (!plotter || !exactSolution) return;

        // 1) серия по h — три уменьшения в 10 раз
        std::vector<double> hs, errs_h;
        double h = h0, tau = tau0;
        for (int k = 0; k < 5; ++k) {
            int N = static_cast<int>(std::round(l / h));
            int M = static_cast<int>(std::round(tMax / tau));

            double err = computeError(kappa, l, tMax, N, M,
                                      f, u0, nu1, nu2, exactSolution);
            hs.push_back(h);
            errs_h.push_back(err);

            h /= 2.0;
        }

        std::vector<std::vector<double>> xs1 = {hs};
        std::vector<std::vector<double>> ys1 = {errs_h};
        std::vector<std::string> labels1 = {"error vs h (sigma=" +
                                            std::to_string(sigma) + ")"};
        plotter->plot(xs1, ys1, labels1, "h", "error");

        // 2) серия по σ при фиксированных шагах
        std::vector<double> sigmas, errs_s;
        for (double s = sigma; s <= 1.0001; s += 0.1) {
            double oldSigma = sigma;
            sigma = s;
            int N = static_cast<int>(std::round(l / h0));
            int M = static_cast<int>(std::round(tMax / tau0));
            double err = computeError(kappa, l, tMax, N, M,
                                      f, u0, nu1, nu2, exactSolution);
            sigma = oldSigma;

            sigmas.push_back(s);
            errs_s.push_back(err);
        }

        std::vector<std::vector<double>> xs2 = {sigmas};
        std::vector<std::vector<double>> ys2 = {errs_s};
        std::vector<std::string> labels2 = {"error vs sigma"};
        plotter->plot(xs2, ys2, labels2, "sigma", "error");
    }

private:
    double computeError(
        double kappa,
        double l, double tMax,
        int N, int M,
        FuncXT f,
        FuncX u0,
        FuncT nu1, FuncT nu2,
        FuncXT exactSolution
    ) {
        const size_t nx  = static_cast<size_t>(N) + 1;
        const size_t nt  = static_cast<size_t>(M) + 1;
        const double h   = l / N;
        const double tau = tMax / M;
        const double mu  = kappa * tau / (h * h);

        Eigen::VectorXd x(nx), t(nt);
        for (size_t i = 0; i < nx; ++i)  x[i] = i * h;
        for (size_t s = 0; s < nt; ++s)  t[s] = s * tau;

        Eigen::MatrixXd u(nt, nx);
        for (size_t i = 0; i < nx; ++i)
            u(0, i) = u0(x[i]);

        for (size_t s = 0; s < nt - 1; ++s) {
            double ts   = t[s];
            double tsp1 = t[s + 1];

            u(s, 0)      = nu1(ts);
            u(s, nx - 1) = nu2(ts);
            u(s + 1, 0)      = nu1(tsp1);
            u(s + 1, nx - 1) = nu2(tsp1);

            if (sigma == 0.0) {
                for (size_t i = 1; i < nx - 1; ++i) {
                    double ui   = u(s, i);
                    double uim1 = u(s, i - 1);
                    double uip1 = u(s, i + 1);
                    double lap  = uip1 - 2.0 * ui + uim1;
                    u(s + 1, i) = ui + mu * lap + tau * f(x[i], ts);
                }
            } else {
                const size_t dim = nx - 2;
                Eigen::MatrixXd A(dim, dim);
                Eigen::VectorXd B(dim);
                A.setZero(); B.setZero();

                for (size_t j = 0; j < dim; ++j) {
                    size_t i = j + 1;

                    if (j > 0)       A(j, j - 1) = -sigma * mu;
                    A(j, j) = 1.0 + 2.0 * sigma * mu;
                    if (j + 1 < dim) A(j, j + 1) = -sigma * mu;

                    double uim1 = u(s, i - 1);
                    double ui   = u(s, i);
                    double uip1 = u(s, i + 1);
                    double lap_n = uip1 - 2.0 * ui + uim1;

                    double rhs = ui + (1.0 - sigma) * mu * lap_n
                                 + tau * f(x[i], ts);

                    if (i == 1)
                        rhs += sigma * mu * u(s + 1, 0);
                    if (i == nx - 2)
                        rhs += sigma * mu * u(s + 1, nx - 1);

                    B(j) = rhs;
                }

                auto result = solver.solve(A, B);
                for (size_t j = 0; j < dim; ++j)
                    u(s + 1, j + 1) = result.solution[j];
            }
        }

        // максимальная ошибка на последнем слое
        double t_last = t[nt - 1];
        double err = 0.0;
        for (size_t i = 0; i < nx; ++i) {
            double diff = std::abs(u(nt - 1, i) - exactSolution(x[i], t_last));
            if (diff > err) err = diff;
        }
        return err;
    }

private:
    ISolver& solver;
    Plotter* plotter;
    double sigma;
    PlotMode mode;
};


#endif // NUMERICAL_METHODS_IN_PHYSICS_TASKHEATEQ_H
