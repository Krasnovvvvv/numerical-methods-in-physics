#ifndef NUMERICAL_METHODS_IN_PHYSICS_INTEGRALEQUATIONTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_INTEGRALEQUATIONTASK_H

#pragma once

#include "Base/IIntegralSolver.h"
#include "Helpers/Timer.h"
#include "Helpers/Plotter.h"

#include <functional>
#include <vector>
#include <string>
#include <optional>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <Eigen/Dense>

class IntegralEquationTask {
public:
    using Kernel   = std::function<double(double,double)>; // Q(t,s)
    using RHS      = std::function<double(double)>; // f(t)
    using Solution = std::function<double(double)>; // exact x(t)

    // mode = 1: график решений; mode = 2: график ошибки |x_num - x_exact|
    explicit IntegralEquationTask(
        IIntegralSolver& integrator,
        Plotter* plotter = nullptr,
        unsigned short mode = 1)
        : integrator(integrator)
        , plotter(plotter)
        , mode(mode)
    {}

    // ---------- Fredholm: x(t) + λ ∫_a^b Q(t,s)x(s)ds = f(t) ----------
    std::vector<double> solve_fredholm(
            double lambda,
            double a, double b,
            std::size_t nodes,
            const Kernel& kernel,
            const RHS& rhs,
            const Solution& exact)
    {
        build_grid(a, b, nodes);
        const std::size_t n = nodes;

        Eigen::MatrixXd A(n, n);
        Eigen::VectorXd b_vec(n);

        auto basis = [&](std::size_t j, double s) {
            static thread_local std::vector<double> e;
            e.assign(n, 0.0);
            e[j] = 1.0;
            return interpolate(s, e);
        };

        for (std::size_t i = 0; i < n; ++i) {
            double ti = grid[i];
            b_vec(static_cast<Eigen::Index>(i)) = rhs(ti);

            for (std::size_t j = 0; j < n; ++j) {
                auto integrand = [&](double s) {
                    return kernel(ti, s) * basis(j, s);
                };
                auto   res= integrator.integrate(integrand, a, b, 0.0);
                double K_ij = res ? res->integral : 0.0;

                A(static_cast<Eigen::Index>(i),
                  static_cast<Eigen::Index>(j)) =
                        (i == j ? 1.0 : 0.0) + lambda * K_ij;
            }
        }

        Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
        Eigen::VectorXd x_vec = lu.solve(b_vec);

        std::vector<double> x(n);
        for (std::size_t i = 0; i < n; ++i)
            x[i] = x_vec(static_cast<Eigen::Index>(i));

        plot_solution("Fredholm: numerical by " + integrator.name() + " (n = " +
            std::to_string(nodes) + ")",
            "Fredholm: exact", grid, x, exact);
        return x;
    }

    // ---------- Volterra: x(t) = rhs(t) + factor ∫_a^t Q(t,s)x(s)ds ----------
    std::vector<double> solve_volterra(
        double a, double b,
        std::size_t   nodes,
        const Kernel& kernel,
        const RHS&    rhs,
        const Solution& exact = Solution(),
        double factor = 1.0)
    {
        build_grid(a, b, nodes);

        std::cout << "\n=== Volterra equation, solver: "
                  << integrator.name() << " ===\n";

        Timer timer;
        std::vector<double> x(nodes, 0.0);

        for (std::size_t i = 0; i < nodes; ++i) {
            double ti = grid[i];

            auto integrand = [&](double s) {
                double xs = interpolate(s, x);
                return kernel(ti, s) * xs;
            };

            auto res= integrator.integrate(integrand, a, ti, 0.0);
            double integral_val = res ? res->integral : 0.0;

            x[i] = rhs(ti) + factor * integral_val;
        }

        auto elapsed_us = timer.elapsed();
        print_solution_stats(grid, x, exact, elapsed_us);
        plot_solution("Volterra: numerical by " + integrator.name() + " (n = " +
            std::to_string(nodes) + ")",
            "Volterra: exact", grid, x, exact);
        return x;
    }

private:
    IIntegralSolver& integrator;
    Plotter* plotter;
    unsigned short mode;
    std::vector<double> grid;

    void build_grid(double a, double b, std::size_t nodes) {
        grid.resize(nodes);
        double h = (b - a) / static_cast<double>(nodes - 1);
        for (std::size_t i = 0; i < nodes; ++i)
            grid[i] = a + h * static_cast<double>(i);
    }

    [[nodiscard]] double interpolate(double s, const std::vector<double>& values) const {
        if (grid.empty())
            return 0.0;
        if (s <= grid.front())
            return values.front();
        if (s >= grid.back())
            return values.back();

        double h = grid[1] - grid[0];
        double t = (s - grid.front()) / h;
        std::size_t k = static_cast<std::size_t>(std::floor(t));
        if (k >= grid.size() - 1)
            return values.back();

        double alpha = t - static_cast<double>(k);
        return (1.0 - alpha) * values[k] + alpha * values[k + 1];
    }

    void print_solution_stats(
        const std::vector<double>& t,
        const std::vector<double>& x,
        const Solution& exact,
        long long elapsed_us) const
    {
        std::cout << "Nodes: " << t.size()
                  << "\nTime: " << elapsed_us << " microseconds\n";

        if (!exact)
            return;

        double max_err = 0.0;
        for (std::size_t i = 0; i < t.size(); ++i)
            max_err = std::max(max_err,
                               std::abs(x[i] - exact(t[i])));

        std::cout << "Max |x_num - x_exact| = "
                  << std::fixed << std::setprecision(8)
                  << max_err << "\n";
    }

    void plot_solution(
        const std::string& num_label,
        const std::string& exact_label,
        const std::vector<double>& t,
        const std::vector<double>& x,
        const Solution& exact) const
    {
        if (!plotter)
            return;

        if (mode == 1) {
            std::vector<std::vector<double>> xs;
            std::vector<std::vector<double>> ys;
            std::vector<std::string> labels;

            xs.push_back(t);
            ys.push_back(x);
            labels.push_back(num_label);

            if (exact) {
                std::vector<double> y_exact(t.size());
                for (std::size_t i = 0; i < t.size(); ++i)
                    y_exact[i] = exact(t[i]);
                xs.push_back(t);
                ys.push_back(std::move(y_exact));
                labels.push_back(exact_label);
            }

            plotter->plot(xs, ys, labels, "t", "x(t)");
        } else if (mode == 2 && exact) {
            std::vector<double> err(t.size());
            for (std::size_t i = 0; i < t.size(); ++i)
                err[i] = std::abs(x[i] - exact(t[i]));

            plotter->plot(t, err,
            "|x_{num} - x_{exact}| for " + integrator.name(),
            "t", "|error|", false);
        }
    }
};

#endif // NUMERICALMETHODS_IN_PHYSICS_INTEGRALEQUATIONTASK_H
