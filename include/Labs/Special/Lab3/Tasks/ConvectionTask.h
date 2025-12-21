#ifndef NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONTASK_H

#pragma once

#include "Labs/Special/Lab3/Base/IConvectionSolver.h"
#include "Helpers/Timer.h"
#include "Helpers/Plotter.h"

#include <vector>
#include <string>
#include <functional>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace special {

struct ConvectionTaskParams {

    // область по x и дискретизация
    double xL;
    double xR;
    std::size_t Nx;

    // область по t и шаг
    double t0;
    double tN;
    double dt;

    // скорость переноса u
    double convectionVelocity;

    // начальное и точное решения
    std::function<double(double)>        u0;    // U(x, t0)
    std::function<double(double,double)> exact; // U_exact(x, t)

    // граничные условия (Дирихле)
    std::function<double(double)> leftBC;  // U(xL, t)
    std::function<double(double)> rightBC; // U(xR, t)
};

class ConvectionTask {

public:
    ConvectionTask(
        IConvectionSolver& solver,
        Plotter* plotter = nullptr,
        std::string x_name = "x",
        std::string u_name = "U"
    )
        : solver(solver),
          plotter(plotter),
          xName(std::move(x_name)),
          uName(std::move(u_name))
    {}

    /// output_times — список моментов, в которые нужно сохранить и нарисовать U(x,t)
    void run(const ConvectionTaskParams& p,
             std::vector<double> output_times)
    {
        // --- дискретизация по x ---
        const std::size_t N  = p.Nx;
        const double dx      = (p.xR - p.xL) / (N - 1);
        const double c       = p.convectionVelocity * p.dt / dx; // c = u dt / dx

        std::sort(output_times.begin(), output_times.end());

        std::vector<double> x(N);
        for (std::size_t i = 0; i < N; ++i)
            x[i] = p.xL + dx * static_cast<double>(i);

        // --- начальное условие ---
        std::vector<double> u_cur(N), u_next(N);
        for (std::size_t i = 0; i < N; ++i)
            u_cur[i] = p.u0 ? p.u0(x[i]) : 0.0;

        // граничные значения на t0
        if (p.leftBC)  u_cur.front() = p.leftBC(p.t0);
        if (p.rightBC) u_cur.back()  = p.rightBC(p.t0);

        // --- подготовка хранения результатов ---
        Timer timer;
        double t = p.t0;
        std::size_t stepCount = 0;

        std::size_t out_idx   = 0;
        const std::size_t out_count = output_times.size();

        std::vector<double>              saved_times;
        std::vector<std::vector<double>> saved_numeric;
        std::vector<std::vector<double>> saved_exact;

        auto save_if_needed = [&](double time, const std::vector<double>& u_vec) {
            while (out_idx < out_count && output_times[out_idx] <= time + 1e-12) {
                double tt = output_times[out_idx];
                saved_times.push_back(tt);

                // численное решение (берём u_vec как ближайшее по времени)
                saved_numeric.push_back(u_vec);

                // точное решение
                std::vector<double> u_ex(N, NAN);
                if (p.exact) {
                    for (std::size_t i = 0; i < N; ++i)
                        u_ex[i] = p.exact(x[i], tt);     // U_exact(x_i, t)
                }
                saved_exact.push_back(u_ex);

                ++out_idx;
            }
        };

        save_if_needed(t, u_cur);

        // --- цикл по времени ---
        while (t < p.tN - 1e-12) {
            double t_next = t + p.dt;
            if (t_next > p.tN) t_next = p.tN;

            // граничные условия на текущем слое
            if (p.leftBC)  u_cur.front() = p.leftBC(t);
            if (p.rightBC) u_cur.back()  = p.rightBC(t);

            solver.step(u_next, u_cur, c);

            // граничные условия на новом слое
            if (p.leftBC)  u_next.front() = p.leftBC(t_next);
            if (p.rightBC) u_next.back()  = p.rightBC(t_next);

            u_cur.swap(u_next);
            t = t_next;
            ++stepCount;

            save_if_needed(t, u_cur);
        }

        auto elapsed_us = timer.elapsed();

        std::cout << "\n=== " << solver.name() << " ===\n";
        std::cout << "Steps: " << stepCount << "\n";
        std::cout << "Computation time: " << elapsed_us << " microseconds\n";
        std::cout << "Convection number c = " << c << "\n";
        std::cout << std::fixed << std::setprecision(7);

        // max‑ошибка по x для каждого сохранённого времени
        std::vector<double> maxErrors(saved_times.size(), NAN);
        if (p.exact) {
            for (std::size_t k = 0; k < saved_times.size(); ++k) {
                double maxErr = 0.0;
                for (std::size_t i = 0; i < N; ++i) {
                    double diff = std::abs(saved_numeric[k][i] - saved_exact[k][i]);
                    if (diff > maxErr) maxErr = diff;
                }
                maxErrors[k] = maxErr;
            }
        }

        for (std::size_t k = 0; k < saved_times.size(); ++k) {
            std::cout << "\n--- t = " << saved_times[k] << " ---\n";
            std::cout << xName << " | Num. " << uName << " | Exact | Error\n";

            for (std::size_t i = 0; i < N; ++i) {
                double u_num = saved_numeric[k][i];
                double u_ex  = p.exact ? saved_exact[k][i] : NAN;
                double err   = p.exact ? std::abs(u_num - u_ex) : NAN;

                std::cout << std::setw(9) << x[i]   << " | "
                          << std::setw(9) << u_num << " | "
                          << std::setw(9) << u_ex  << " | "
                          << std::setw(9) << err   << "\n";
            }

            if (p.exact)
                std::cout << "Max error at this time: " << maxErrors[k] << "\n";
        }

        // --- графики ---
        if (plotter && !saved_times.empty()) {
            // Один график: все моменты времени + точное решение
            std::vector<std::vector<double>> xs;
            std::vector<std::vector<double>> ys;
            std::vector<std::string>         labels;

            for (std::size_t k = 0; k < saved_times.size(); ++k) {
                // численное решение в момент t_k
                xs.push_back(x);
                ys.push_back(saved_numeric[k]);
                labels.push_back(
                    uName + " numeric t=" + std::to_string(saved_times[k]) +
                    " by " + solver.name()
                );

                if (p.exact) {
                    xs.push_back(x);
                    ys.push_back(saved_exact[k]);
                    labels.push_back(
                        uName + " exact t=" + std::to_string(saved_times[k])
                    );
                }
            }

            plotter->plot(xs, ys, labels,
                          xName,
                          uName + "(x, various t)");

            // график max‑ошибки по времени
            if (p.exact) {
                std::vector<std::vector<double>> xs_err{ saved_times };
                std::vector<std::vector<double>> ys_err{ maxErrors };
                std::vector<std::string>         labels_err{ "Max error over x" };

                plotter->plot(xs_err, ys_err, labels_err,
                              "t",
                              "Max |U_{num} - U_{exact}|");
            }
        }
    }

private:
    IConvectionSolver& solver;
    Plotter*           plotter;
    std::string        xName;
    std::string        uName;
};

} // namespace special

#endif // NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONTASK_H