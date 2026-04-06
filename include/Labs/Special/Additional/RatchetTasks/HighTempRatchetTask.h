#ifndef NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPRATCHETTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPRATCHETTASK_H

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
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
    unsigned int run_seed_;

public:
    HighTempRatchetTask(Plotter& plotter,
                        RatchetParams& params,
                        std::size_t n_particles,
                        std::size_t total_time,
                        bool build_plot = true,
                        std::size_t max_plot_points = 800,
                        unsigned int run_seed = 1u)
        : plotter_(plotter),
          params_(params),
          n_particles_(n_particles),
          total_time_(total_time),
          build_plot_(build_plot),
          max_plot_points_(max_plot_points),
          run_seed_(run_seed) {}

    void run() {
    auto F1 = [](const double xi) -> double {
        return (3.0 * xi * (1.0 + 2.0 * xi)) /
               ((1.0 + 4.0 * xi) * (1.0 + 4.0 * xi) * (1.0 + xi));
    };

    const double xi = 1.0 / (4.0 * M_PI * M_PI * params_.epsilon);

    const double v_analytical =
        0.25 * M_PI * params_.V1 * params_.V1 * params_.V2 *
        (1.0 - params_.alpha) * (1.0 - params_.alpha) * (1.0 + params_.alpha) *
        F1(xi);

    const double v_analytical5 = highTempVelocity5(params_);

    std::cout << std::setprecision(12)
              << "Velocity (analytical): " << v_analytical << "\n";

    std::cout << std::setprecision(12)
              << "Velocity (analytical, refined): " << v_analytical5 << "\n";

    const unsigned int base_seed_dicho = mix_seed_
    (
        600u, run_seed_, 0x13579BDFu
    );
    const unsigned int base_seed_gauss = mix_seed_
    (
        1600u, run_seed_, 0x2468ACE1u
    );

    DimLessLangevinSolver solver(params_, base_seed_gauss);

    std::vector<DichotomicNoise> noises;
    noises.reserve(n_particles_);

    for (std::size_t p = 0; p < n_particles_; ++p) {
        noises.emplace_back(
            params_.a,
            params_.epsilon,
            params_.dt,
            make_particle_seed_(base_seed_dicho, p, 0xA5A5A5A5u)
        );
    }

    const auto N = static_cast<std::size_t>(total_time_ / params_.dt);
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
    const double sol_v = solution.mean_velocity;

    const double error_abs = std::abs(sol_v - v_analytical);
    const double error_rel = (v_analytical != 0.0)
        ? (error_abs / std::abs(v_analytical)) * 100.0
        : 100.0;

    std::cout << std::setprecision(12)
              << "Velocity (numerical, regression): " << sol_v << "\n";

    std::cout << "Abs error: " << error_abs << "\n"
              << "Relative error: " << error_rel << " %\n"
              << "Simulation time: " << elapsed / 60.0 << " min\n";

    if (!solution.has_trajectory
        || solution.sim_time.size() < 2
        || solution.mean_x.size() < 2) {
        return;
    }

    const std::vector<double>& t_plot = solution.sim_time;
    const std::vector<double>& x_plot = solution.mean_x;

    const std::size_t match_start =
        find_match_start_index_(t_plot, x_plot,
                                sol_v,
                                v_analytical5);

    const double t_match = t_plot[match_start];
    const double x_match = x_plot[match_start];

    std::cout << "Match start index: " << match_start << "\n"
              << "Match start time: " << t_match << "\n";

    constexpr double fallback_trend_shift = 0.0;

    double trend_shift = fallback_trend_shift;
    if (std::isfinite(t_match) && std::isfinite(x_match)) {
        trend_shift = x_match - sol_v * t_match;
    }
    if (!std::isfinite(trend_shift)) {
        trend_shift = fallback_trend_shift;
    }

    std::cout << "Trend vertical shift: " << trend_shift << "\n";

    std::vector<double> t_tail(
        t_plot.begin() + static_cast<std::ptrdiff_t>(match_start),
        t_plot.end()
    );

    std::vector<double> x_an;
    std::vector<double> x_an_refined;
    std::vector<double> x_trend;

    x_an.reserve(t_plot.size());
    x_an_refined.reserve(t_plot.size());
    x_trend.reserve(t_tail.size());

    for (double t : t_plot) {
        x_an.push_back(v_analytical * t);
        x_an_refined.push_back(v_analytical5 * t);
    }

    for (double t : t_tail) {
        x_trend.push_back(trend_shift + sol_v * t);
    }

    std::vector<std::vector<double>> xs;
    std::vector<std::vector<double>> ys;
    std::vector<std::string> labels;

    xs.push_back(t_plot);
    ys.push_back(x_plot);
    labels.push_back("Numerical");

    xs.push_back(t_plot);
    ys.push_back(x_an);
    labels.push_back("Analytical");

    xs.push_back(t_plot);
    ys.push_back(x_an_refined);
    labels.push_back("Analytical (refined)");

    xs.push_back(t_tail);
    ys.push_back(x_trend);
    labels.push_back("Trend");

    plotter_.plot(xs, ys, labels, "t/{/Symbol t}_{D}", "{x / L}",
                  {false, false, false, true});
}

private:
    static std::size_t find_match_start_index_(const std::vector<double>& t,
                                               const std::vector<double>& x,
                                               double v_numerical,
                                               double v_refined)
    {
        const std::size_t n = std::min(t.size(), x.size());
        if (n < 8) {
            return 0;
        }

        const std::size_t half_window = std::max<std::size_t>(2, n / 40);
        const std::size_t confirm_count = std::max<std::size_t>(3, n / 30);
        const std::size_t start_search = std::min<std::size_t>(n - 1, n / 10);

        auto local_slope = [&](std::size_t i) -> double {
            const std::size_t left = (i > half_window) ? (i - half_window) : 0;
            const std::size_t right = std::min<std::size_t>(n - 1, i + half_window);

            const double dt = t[right] - t[left];
            if (dt <= 0.0) {
                return std::numeric_limits<double>::quiet_NaN();
            }

            return (x[right] - x[left]) / dt;
        };

        const double abs_tol_num = std::max(1e-6, 0.10 * std::abs(v_numerical));
        const double abs_tol_ref = std::max(1e-6, 0.10 * std::abs(v_refined));

        std::size_t streak = 0;

        for (std::size_t i = start_search + half_window; i + half_window < n; ++i) {
            const double slope = local_slope(i);
            if (!std::isfinite(slope)) {
                streak = 0;
                continue;
            }

            const bool close_to_num = std::abs(slope - v_numerical) <= abs_tol_num;
            const bool close_to_ref = std::abs(slope - v_refined) <= abs_tol_ref;

            if (close_to_num && close_to_ref) {
                ++streak;
                if (streak >= confirm_count) {
                    return i - streak + 1;
                }
            } else {
                streak = 0;
            }
        }

        return start_search;
    }

    static unsigned int mix_seed_(unsigned int base,
                                  unsigned int run_seed,
                                  std::uint32_t tag)
    {
        std::seed_seq seq{
            static_cast<std::uint32_t>(base),
            static_cast<std::uint32_t>(run_seed),
            static_cast<std::uint32_t>(tag),
            0x9E3779B9u,
            0x85EBCA6Bu,
            0xC2B2AE35u
        };

        std::uint32_t out[1];
        seq.generate(std::begin(out), std::end(out));
        return static_cast<unsigned int>(out[0]);
    }

    static unsigned int make_particle_seed_(unsigned int base_seed,
                                            std::size_t particle_index,
                                            std::uint32_t tag)
    {
        const auto p64 = static_cast<std::uint64_t>(particle_index);

        std::seed_seq seq{
            static_cast<std::uint32_t>(base_seed),
            static_cast<std::uint32_t>(tag),
            static_cast<std::uint32_t>(p64 & 0xFFFFFFFFu),
            static_cast<std::uint32_t>((p64 >> 32) & 0xFFFFFFFFu),
            0x517CC1B7u,
            0xD2B74407u
        };

        std::uint32_t out[1];
        seq.generate(std::begin(out), std::end(out));
        return static_cast<unsigned int>(out[0]);
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPRATCHETTASK_H