#ifndef NUMERICAL_METHODS_IN_PHYSICS_DIFFUSION_VISUALIZER_H
#define NUMERICAL_METHODS_IN_PHYSICS_DIFFUSION_VISUALIZER_H

#pragma once

#include "BaseDiffusionSolver.h"
#include "DiffusionParameters.h"
#include "Potential1D.h"
#include "Helpers/Plotter.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace special {

/**
 * @class DiffusionVisualizer
 * @brief Визуализация результатов диффузии через Plotter с matplot++
 *
 * РЕАЛЬНАЯ отрисовка графиков через Plotter::plot()
 */
class DiffusionVisualizer {
public:
    explicit DiffusionVisualizer(Plotter* plotter = nullptr)
        : plotter_(plotter) {}

    // ========================================================================
    // СЕМЕЙСТВО ТРАЕКТОРИЙ (несколько частиц на одном графике)
    // ========================================================================

    void plot_trajectories(const DiffusionResult& result,
                          int max_particles = 10,
                          const std::string& title = "Particle Trajectories") const {
        if (result.trajectories.empty() || !plotter_) return;

        int n_plot = std::min(max_particles, static_cast<int>(result.trajectories.size()));

        std::vector<std::vector<double>> xs_traj;
        std::vector<std::vector<double>> ys_traj;
        std::vector<std::string> labels_traj;

        for (int p = 0; p < n_plot; ++p) {
            xs_traj.push_back(result.t);
            ys_traj.push_back(result.trajectories[p]);
            labels_traj.push_back("Particle " + std::to_string(p + 1));
        }

        plotter_->plot(xs_traj, ys_traj, labels_traj, "Time t", "Position x(t)");
        std::cout << "  ✓ " << title << " (" << n_plot << " particles plotted)\n";
    }

    // ========================================================================
    // ВСЕ 4 МОМЕНТА: m₁, m₂, m₃, m₄
    // ========================================================================

    void plot_all_moments(const DiffusionResult& result,
                         const DiffusionParameters& params) const {
        if (result.t.empty() || !plotter_) return;

        // ГРАФИК 1: m₁ = <x>
        if (!result.mean_x.empty()) {
            plotter_->plot(result.t, result.mean_x,
                "m₁ = <x>", "Time t", "First Moment <x>");
            std::cout << "  ✓ First moment m₁ = <x>\n";
        }

        // ГРАФИК 2: m₂ = <x²> vs теория
        if (!result.mean_x2.empty()) {
            std::vector<double> theory_m2(result.t.size());
            for (size_t i = 0; i < result.t.size(); ++i) {
                theory_m2[i] = 2.0 * params.diffusion_coeff * result.t[i];
            }

            std::vector<std::vector<double>> xs_m2 = {result.t, result.t};
            std::vector<std::vector<double>> ys_m2 = {result.mean_x2, theory_m2};
            std::vector<std::string> labels_m2 = {"Numerical <x²>", "Theory 2Dt"};
            plotter_->plot(xs_m2, ys_m2, labels_m2, "Time t", "Second Moment <x²>");
            std::cout << "  ✓ Second moment m₂ = <x²> vs theory\n";
        }

        // ГРАФИК 3: m₃ = <x³>
        if (!result.mean_x3.empty()) {
            plotter_->plot(result.t, result.mean_x3,
                "m₃ = <x³>", "Time t", "Third Moment <x³>");
            std::cout << "  ✓ Third moment m₃ = <x³>\n";
        }

        // ГРАФИК 4: m₄ = <x⁴>
        if (!result.mean_x4.empty()) {
            plotter_->plot(result.t, result.mean_x4,
                "m₄ = <x⁴>", "Time t", "Fourth Moment <x⁴>");
            std::cout << "  ✓ Fourth moment m₄ = <x⁴>\n";
        }
    }

    // ========================================================================
    // ГИСТОГРАММА (через Plotter если поддерживает)
    // ========================================================================

    void plot_distribution(const DiffusionResult& result,
                          int n_bins = 50,
                          const std::string& title = "Final Distribution") const {
        if (result.trajectories.empty() || !plotter_) return;

        std::vector<double> final_positions;
        for (const auto& traj : result.trajectories) {
            if (!traj.empty()) {
                final_positions.push_back(traj.back());
            }
        }

        if (final_positions.empty()) return;

        // Сортируем позиции для гистограммы
        std::sort(final_positions.begin(), final_positions.end());

        // Вычисляем бины
        double min_x = final_positions.front();
        double max_x = final_positions.back();
        double bin_width = (max_x - min_x) / n_bins;
        if (bin_width < 1e-10) bin_width = 1.0;

        std::vector<double> bin_centers(n_bins);
        std::vector<double> bin_counts(n_bins, 0);

        for (double pos : final_positions) {
            int bin = static_cast<int>((pos - min_x) / bin_width);
            if (bin >= n_bins) bin = n_bins - 1;
            if (bin < 0) bin = 0;
            bin_counts[bin]++;
        }

        for (int i = 0; i < n_bins; ++i) {
            bin_centers[i] = min_x + (i + 0.5) * bin_width;
        }

        plotter_->plot(bin_centers, bin_counts,
            "Distribution", "Position x", "Frequency");
        std::cout << "  ✓ " << title << " (histogram)\n";
    }

    // ========================================================================
    // СРАВНЕНИЕ МСД С ТЕОРИЕЙ
    // ========================================================================

    void compare_with_theory(const DiffusionResult& result,
                            const DiffusionParameters& params) const {
        if (result.t.empty() || !plotter_) return;

        std::vector<double> theory(result.t.size());
        for (size_t i = 0; i < result.t.size(); ++i) {
            theory[i] = 2.0 * params.diffusion_coeff * result.t[i];
        }

        std::vector<std::vector<double>> xs_msd = {result.t, result.t};
        std::vector<std::vector<double>> ys_msd = {result.displacement_squared, theory};
        std::vector<std::string> labels_msd = {"Numerical MSD", "Theory 2Dt"};

        plotter_->plot(xs_msd, ys_msd, labels_msd, "Time t", "<(Δx)²>");
        std::cout << "  ✓ MSD comparison with theory\n";
    }

    // ========================================================================
    // СРЕДНЕЕ ПОЛОЖЕНИЕ
    // ========================================================================

    void plot_mean_position(const DiffusionResult& result,
                           const std::string& title = "Mean Position") const {
        if (result.t.empty() || !plotter_) return;

        plotter_->plot(result.t, result.mean_x,
            title, "Time t", "<x(t)>");
        std::cout << "  ✓ " << title << "\n";
    }

    // ========================================================================
    // ПОТЕНЦИАЛ U(x)
    // ========================================================================

    void plot_potential(const Potential1D& potential,
                       const DiffusionParameters& params,
                       const std::string& title = "Potential U(x)") const {
        if (!plotter_) return;

        std::vector<double> x;
        std::vector<double> U;
        int n_points = 500;

        for (int i = 0; i <= n_points; ++i) {
            double xi = params.x_min + (params.x_max - params.x_min) * i / n_points;
            x.push_back(xi);
            U.push_back(potential.U(xi));
        }

        plotter_->plot(x, U, title, "Position x", "Potential U(x)");
        std::cout << "  ✓ " << title << "\n";
    }

    // ========================================================================
    // ДОПОЛНИТЕЛЬНЫЕ ГРАФИКИ
    // ========================================================================

    void plot_statistics(const DiffusionResult& result,
                        const std::string& main_title = "Diffusion Statistics") const {
        if (!plotter_) return;

        std::cout << "  ✓ " << main_title << "\n";

        // Выводим основные показатели
        if (!result.t.empty() && !result.mean_x.empty()) {
            plotter_->plot(result.t, result.mean_x,
                "Mean position", "Time t", "<x>");
        }

        if (!result.t.empty() && !result.displacement_squared.empty()) {
            plotter_->plot(result.t, result.displacement_squared,
                "MSD", "Time t", "<(x-x0)²>");
        }
    }

    void plot_higher_moments(const DiffusionResult& result) const {
        if (!plotter_) return;

        if (!result.t.empty() && !result.mean_x3.empty()) {
            plotter_->plot(result.t, result.mean_x3,
                "m₃ = <x³>", "Time t", "Third Moment");
            std::cout << "  ✓ Third moment\n";
        }

        if (!result.t.empty() && !result.mean_x4.empty()) {
            plotter_->plot(result.t, result.mean_x4,
                "m₄ = <x⁴>", "Time t", "Fourth Moment");
            std::cout << "  ✓ Fourth moment\n";
        }
    }

private:
    Plotter* plotter_;
};

} // namespace special

#endif // NUMERICAL_METHODS_IN_PHYSICS_DIFFUSION_VISUALIZER_H