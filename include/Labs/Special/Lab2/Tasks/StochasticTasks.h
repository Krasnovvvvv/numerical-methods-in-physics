#ifndef NUMERICAL_METHODS_IN_PHYSICS_STOCHASTIC_TASKS_H
#define NUMERICAL_METHODS_IN_PHYSICS_STOCHASTIC_TASKS_H

#pragma once

#include "Labs/Special/Lab2/Base/StochasticTask.h"
#include "Labs/Special/Lab2/Base/Potential1D.h"
#include <vector>
#include <memory>
#include <iostream>

namespace special {

// ============================================================================
// ЗАДАЧА 1 и 3: СВОБОДНАЯ ДИФФУЗИЯ
// ============================================================================

    class FreeDiffusion : public StochasticTask {
    public:
        FreeDiffusion(BaseDiffusionSolver& solver,
                      Plotter* plotter,
                      unsigned short task_id,
                      const std::string& task_name)
            : StochasticTask(solver, plotter, task_id, task_name) {}

    protected:
        std::unique_ptr<Potential1D> create_potential(
                const DiffusionParameters& /*params*/) override {
            // F(x,t) = 0  → свободная диффузия
            return std::make_unique<FreePotential>();
        }

        void visualize(const DiffusionResult& result,
                       const DiffusionParameters& params) override {
            // 1) Демонстрация роста неопределенности: гистограммы + теория p(x,t)
            visualizer_.plot_time_sliced_histograms_with_theory(
                result, params, 4, 50);

            // 2) Семейство траекторий нескольких частиц
            visualizer_.plot_trajectories(
                result,
                10,
                "Trajectories of particles (free diffusion)");

            // 3) МСД численное vs теория 2 D t
            visualizer_.compare_with_theory(result, params);

            // 4) График <x(t)>, <[x(t)-x(0)]^2> и оценка D из наклона
            visualizer_.plot_statistics(result, "Free diffusion: statistics");
        }
    };

// ============================================================================
// ЗАДАЧА 2: АНАЛИЗ МОМЕНТОВ
// ============================================================================

    class MomentsAnalysis : public StochasticTask {
    public:
        MomentsAnalysis(BaseDiffusionSolver& solver,
                        Plotter* plotter,
                        unsigned short task_id,
                        const std::string& task_name)
            : StochasticTask(solver, plotter, task_id, task_name) {}

    protected:
        std::unique_ptr<Potential1D> create_potential(
                const DiffusionParameters& /*params*/) override {
            // Свободная диффузия: F(x,t) = 0
            return std::make_unique<FreePotential>();
        }

        void visualize(const DiffusionResult& result,
                       const DiffusionParameters& params) override {
            // Графики m1(t), m2(t), m3(t), m4(t) с теоретическими кривыми
            visualizer_.plot_all_moments(result, params);

            // оценка того, насколько m1 и m3 близки к нулю при данном числе частиц.
            visualizer_.plot_statistics(result,
                "Moments analysis: check <x> and <x^3> ~ 0");
        }
    };

// ============================================================================
// ЗАДАЧА 4: СТАЦИОНАРНАЯ ОДНОРОДНАЯ СИЛА
// ============================================================================

    class ConstantForceDiffusion : public StochasticTask {
    public:
        ConstantForceDiffusion(BaseDiffusionSolver& solver,
                               Plotter* plotter,
                               unsigned short task_id,
                               const std::string& task_name)
            : StochasticTask(solver, plotter, task_id, task_name) {}

    protected:
        std::unique_ptr<Potential1D> create_potential(
                const DiffusionParameters& params) override {
            return std::make_unique<ConstantForcePotential>(params.constant_force);
        }

        void visualize(const DiffusionResult& result,
                       const DiffusionParameters& params) override {
            ConstantForcePotential pot(params.constant_force);

            visualizer_.plot_potential(
                pot, params, "U(x) = -F x (constant force)");

            visualizer_.plot_trajectories(
                result, 10, "Trajectories with constant force");

            visualizer_.plot_mean_position(
                result, "Mean position under constant force");

            visualizer_.plot_statistics(
                result, "Statistics with constant force");
        }
    };

    // ============================================================================
    // ЗАДАЧА 5A: ПЕРИОДИЧЕСКИЙ ПОТЕНЦИАЛ
    // U(x) = V(x) = V0 sin(2π x / L)
    // ============================================================================

    class PeriodicPotentialDiffusion : public StochasticTask {
    public:
        PeriodicPotentialDiffusion(BaseDiffusionSolver& solver,
                                   Plotter* plotter,
                                   unsigned short task_id,
                                   const std::string& task_name,
                                   double V0,
                                   double L)
            : StochasticTask(solver, plotter, task_id, task_name)
            , V0_(V0)
            , L_(L) {}

    protected:
        std::unique_ptr<Potential1D> create_potential(
                const DiffusionParameters& /*params*/) override {
            // Амплитуда задаётся при создании задачи
            return std::make_unique<SinPotential>(V0_, L_);
        }

        void visualize(const DiffusionResult& result,
                       const DiffusionParameters& params) override {
            SinPotential pot(V0_, L_);

            visualizer_.plot_potential(
                pot, params, "U(x) = V0 * sin²(πx/L)");

            visualizer_.plot_trajectories(
                result, 10, "Trajectories in periodic potential");

            visualizer_.plot_distribution(
                result, 40, "Final distribution in periodic potential");

            visualizer_.plot_statistics(
                result, "Statistics in periodic potential");
        }

    private:
        double V0_;
        double L_;
    };

    // ============================================================================
    // ЗАДАЧА 5B: ПЕРИОДИЧЕСКИЙ ПОТЕНЦИАЛ + ПОСТОЯННАЯ СИЛА
    // U(x,t) = V(x) - F x,  V(x) = V0 sin(2π x / L)
    // ============================================================================

    class PeriodicPotentialWithForce : public StochasticTask {
    public:
        PeriodicPotentialWithForce(BaseDiffusionSolver& solver,
                                   Plotter* plotter,
                                   unsigned short task_id,
                                   const std::string& task_name,
                                   double V0,
                                   double L,
                                   double F_static)
            : StochasticTask(solver, plotter, task_id, task_name)
            , V0_(V0)
            , L_(L)
            , F_static_(F_static) {}

    protected:
        std::unique_ptr<Potential1D> create_potential(
                const DiffusionParameters& /*params*/) override {
            auto base = std::make_unique<SinPotential>(V0_, L_);
            // F_amplitude = 0 → только стационарный наклон
            return std::make_unique<TiltingPotential>(
                std::move(base), F_static_, /*F_amplitude=*/0.0, /*omega=*/1.0);
        }

        void visualize(const DiffusionResult& result,
                       const DiffusionParameters& params) override {
            auto base = std::make_unique<SinPotential>(V0_, L_);
            TiltingPotential pot(std::move(base), F_static_, 0.0, 1.0);

            visualizer_.plot_potential(
                pot, params, "U(x) = V0 sin(2π x / L) - F x");

            visualizer_.plot_trajectories(
                result, 10, "Trajectories in periodic potential with force");

            visualizer_.plot_mean_position(
                result, "Mean position (directed transport)");

            visualizer_.plot_statistics(
                result, "Statistics in periodic potential + constant force");
        }

    private:
        double V0_;
        double L_;
        double F_static_;
    };

} // namespace special

#endif // NUMERICAL_METHODS_IN_PHYSICS_STOCHASTIC_TASKS_H