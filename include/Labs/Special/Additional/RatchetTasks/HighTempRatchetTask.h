#ifndef NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPRATCHETTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPRATCHETTASK_H

#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Labs/Special/Additional/DffusionSolvers/DimLessLangevinSolver.h"
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"
#include "Labs/Special/Additional/Base/RatchetParams.h"
#include "Labs/Special/Additional/Helpers/HighTempVelocity5.h"
#include "Helpers/Plotter.h"
#include "Helpers/Timer.h"

class HighTempRatchetTask {
    Plotter& plotter_;
    RatchetParams& params_;
    std::size_t n_particles_;
    std::size_t total_time_;
    bool build_plot_;
    std::size_t max_plot_points_;

public:
    HighTempRatchetTask(Plotter& plotter,
                        RatchetParams& params,
                        std::size_t n_particles,
                        std::size_t total_time,
                        bool build_plot = true,
                        std::size_t max_plot_points = 800)
        : plotter_(plotter),
          params_(params),
          n_particles_(n_particles),
          total_time_(total_time),
          build_plot_(build_plot),
          max_plot_points_(max_plot_points) {}

    void run() {
        auto F1 = [](const double xi) -> double {
            return (3.0 * xi * (1.0 + 2.0 * xi)) /
                   ((1.0 + 4.0 * xi) * (1.0 + 4.0 * xi) * (1.0 + xi));
        };

        const double xi = 1.0 / (4.0 * M_PI * M_PI * params_.epsilon);
        double v_analytical =
            0.25 * M_PI * params_.V1 * params_.V1 * params_.V2 *
            (1.0 - params_.alpha) * (1.0 - params_.alpha) * (1.0 + params_.alpha) *
            F1(xi);

        const double v_analytical5 = highTempVelocity5(params_);

        std::cout << std::setprecision(12)
                  << "Velocity (analytical): " << v_analytical << "\n";

        std::cout << std::setprecision(12)
          << "Velocity (analytical, adjusted): " << v_analytical5 << "\n";

        /*if (v_analytical5 < v_analytical)
            v_analytical = v_analytical5;*/

        DimLessLangevinSolver solver(params_, 600u);

        std::vector<DichotomicNoise> noises;
        noises.reserve(n_particles_);
        for (std::size_t p = 0; p < n_particles_; ++p) {
            noises.emplace_back(params_.a, params_.epsilon, params_.dt,
                                600u + static_cast<unsigned int>(p));
        }

        const std::size_t N = static_cast<std::size_t>(total_time_ / params_.dt);
        const std::size_t burn_in = (N > 0) ? (N / 10) : 0;

        const bool store_trajectory = build_plot_;
        const std::size_t trajectory_stride =
            store_trajectory
                ? std::max<std::size_t>(1, N / std::max<std::size_t>(1, max_plot_points_))
                : 1;

        Timer<std::chrono::duration<double>> timer;
        auto solution = solver.solve(
            noises,
            n_particles_,
            total_time_,
            burn_in,
            0.0,
            store_trajectory,
            trajectory_stride
        );
        const double elapsed = timer.elapsed();

        double v_full = 0.0;
        if (solution.has_trajectory &&
            solution.sim_time.size() >= 2 &&
            solution.mean_x.size() == solution.sim_time.size())
        {
            const double dt_full = solution.sim_time.back() - solution.sim_time.front();
            if (dt_full > 0.0) {
                v_full = (solution.mean_x.back() - solution.mean_x.front()) / dt_full;
            }
        }

        const double error_abs = std::abs(solution.mean_velocity - v_analytical);
        const double error_rel = (v_analytical != 0.0)
            ? (error_abs / std::abs(v_analytical)) * 100.0
            : 100.0;

        std::cout << std::setprecision(12)
                  << "Velocity (numerical, regression): " << solution.mean_velocity << "\n";

        if (solution.has_trajectory) {
            std::cout << "Velocity (full-window): " << v_full << "\n";
        }

        std::cout << "Abs error: " << error_abs << "\n"
                  << "Relative error: " << error_rel << " %\n"
                  << "Simulation time: " << elapsed << " sec (" << elapsed / 60.0 << " min)\n";

        if (!solution.has_trajectory || solution.sim_time.size() < 2 || solution.mean_x.size() < 2) {
            return;
        }

        std::vector<double> t_plot = solution.sim_time;
        std::vector<double> x_plot = solution.mean_x;
        std::vector<double> x_an, x_an_corr;
        x_an.reserve(t_plot.size());
        x_an_corr.reserve(t_plot.size());

        const double t0 = t_plot.front();
        const double x0u = x_plot.front();
        for (double t : t_plot) {
            x_an.push_back(x0u + v_analytical * (t - t0));
            x_an_corr.push_back(x0u + v_analytical5 * (t - t0));
        }

        std::vector<std::vector<double>> xs = {t_plot, t_plot, t_plot};
        std::vector<std::vector<double>> ys = {x_plot, x_an, x_an_corr};
        std::vector<std::string> labels = {
            "Numerical",
            "Analytical",
            "Analytical (adjusted)"
        };

        plotter_.plot(xs, ys, labels, "t/{tau_D}", "{x / L}");
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPRATCHETTASK_H