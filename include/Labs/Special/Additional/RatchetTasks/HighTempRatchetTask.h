#ifndef NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPRATCHETTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPRATCHETTASK_H

#pragma once
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "Labs/Special/Additional/DffusionSolvers/DimLessLangevinSolver.h"
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"
#include "Labs/Special/Additional/Base/RatchetParams.h"
#include "Helpers/Plotter.h"
#include "Helpers/Timer.h"

class HighTempRatchetTask {
    Plotter& plotter_;
    RatchetParams& params_;
    std::size_t n_particles_;
    std::size_t total_time_;

public:
    HighTempRatchetTask(
        Plotter& plotter,
        RatchetParams& params,
        std::size_t n_particles,
        std::size_t total_time)
        :   plotter_(plotter),
            params_(params),
            n_particles_(n_particles),
            total_time_(total_time){}

    void run() {
        auto F1 = [](const double xi) -> double {
            return (3.0 * xi * (1.0 + 2.0 * xi)) /
                ((1 + 4 * xi) * (1 + 4 * xi) * (1 + xi));
        };

        const double xi = 1.0 / (4.0 * M_PI * M_PI * params_.epsilon);
        const double v_analytical = 0.25 * M_PI * params_.V1 * params_.V1 *
            params_.V2 * (1.0 - params_.alpha) * (1.0 - params_.alpha) *
                (1.0 + params_.alpha) * F1(xi);
        std::cout << std::setprecision(12)
                  << "Velocity (analytical): " << v_analytical << "\n";

        DimLessLangevinSolver solver(params_, 600u);

        std::vector<DichotomicNoise> noises;
        noises.reserve(n_particles_);
        for (std::size_t p = 0; p < n_particles_; ++p) {
            noises.emplace_back(params_.a,
                params_.epsilon, params_.dt,
                600u +static_cast<unsigned int>(p));
        }

        Timer<std::chrono::duration<double>> timer;
        auto solution = solver.solve(noises, n_particles_, total_time_);
        const double elapsed = timer.elapsed();

        const auto& mean_x = solution.mean_x;

        double v_full = 0.0;
        if (solution.sim_time.size() >= 2 && mean_x.size() == solution.sim_time.size()) {
            const double dt_full = solution.sim_time.back() - solution.sim_time.front();
            if (dt_full > 0.0) {
                v_full = (mean_x.back() - mean_x.front()) / dt_full;
            }
        }

        const double error_abs = std::abs(v_full - v_analytical);
        const double error_rel = (v_analytical != 0.0)
            ? (error_abs / std::abs(v_analytical)) * 100.0
            : 100.0;

        std::cout << std::setprecision(12)
                  << "Velocity (numerical): " << solution.mean_velocity << "\n"
                  << "Velocity (full-window): " << v_full << "\n"
                  << "Abs error: " << error_abs << "\n"
                  << "Relative error: " << error_rel << " %\n"
                  << "Simulation time: " << elapsed << " sec (" << elapsed / 60.0 << " min)\n";

        const std::size_t N = static_cast<std::size_t>(
                total_time_ / params_.dt);
        const std::size_t step = std::max<std::size_t>(1, N / 800);

        std::vector<double> t_plot, x_plot, x_an;
        t_plot.reserve(N / step + 1);
        x_plot.reserve(N / step + 1);
        x_an.reserve(N / step + 1);

        const std::size_t i0 = std::min<std::size_t>(
        N / 10,
        (solution.sim_time.size() > 0 ? (solution.sim_time.size() - 1) : 0)
        );

        const double t0  = solution.sim_time.empty() ? 0.0 : solution.sim_time[i0];
        const double x0u = mean_x.empty() ? 0.0 : mean_x[i0];

        for (std::size_t i = i0;
            i < solution.sim_time.size() && i < mean_x.size(); i += step) {
            const double t = solution.sim_time[i];
            t_plot.push_back(t);
            x_plot.push_back(mean_x[i]);
            x_an.push_back(x0u + v_analytical * (t - t0));
        }

        std::vector<std::vector<double>> xs = {t_plot, t_plot};
        std::vector<std::vector<double>> ys = {x_plot, x_an};
        std::vector<std::string> labels = {
            "Numerical",
            "Analytical"
        };

        plotter_.plot(xs, ys, labels, "t/{tau_D}", "{x / L}");
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPRATCHETTASK_H
