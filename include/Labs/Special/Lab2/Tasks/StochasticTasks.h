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
// ЗАДАЧА 1: СВОБОДНАЯ ДИФФУЗИЯ
// ============================================================================

class FreeDiffusion : public StochasticTask {
public:
    FreeDiffusion(BaseDiffusionSolver& solver,
                 Plotter* plotter,
                 unsigned short task_id,
                 const std::string& task_name)
        : StochasticTask(solver, plotter, task_id, task_name) {}

protected:
    std::unique_ptr<Potential1D> create_potential(const DiffusionParameters& params) override {
        return std::make_unique<FreePotential>();
    }

    void visualize(const DiffusionResult& result,
                  const DiffusionParameters& params) override {
        visualizer_.plot_time_sliced_histograms_with_theory(result, params, 4);
        visualizer_.plot_trajectories(result, 5,
            "Particle trajectories (free diffusion)");
        visualizer_.compare_with_theory(result, params);
        visualizer_.plot_statistics(result, "Free Diffusion Statistics");
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
    std::unique_ptr<Potential1D> create_potential(const DiffusionParameters& params) override {
        return std::make_unique<FreePotential>();
    }

    void visualize(const DiffusionResult& result,
                  const DiffusionParameters& params) override {
        visualizer_.plot_all_moments(result, params);
    }
};

// ============================================================================
// ЗАДАЧА 3: ДИФФУЗИЯ В ПОТЕНЦИАЛЕ
// ============================================================================

class DiffusionInPotential : public StochasticTask {
public:
    DiffusionInPotential(BaseDiffusionSolver& solver,
                        Plotter* plotter,
                        unsigned short task_id,
                        const std::string& task_name)
        : StochasticTask(solver, plotter, task_id, task_name) {}

protected:
    std::unique_ptr<Potential1D> create_potential(const DiffusionParameters& params) override {
        return std::make_unique<HarmonicPotential>(params.potential_stiffness);
    }

    void visualize(const DiffusionResult& result,
                  const DiffusionParameters& params) override {
        HarmonicPotential pot(params.potential_stiffness);
        visualizer_.plot_potential(pot, params, "Harmonic potential U(x)");
        visualizer_.plot_trajectories(result, 15, "Trajectories (relaxation)");
        visualizer_.plot_distribution(result, 30, "Final distribution (Boltzmann)");
        visualizer_.plot_statistics(result, "Diffusion in Potential");
    }
};

// ============================================================================
// ЗАДАЧА 4: МИГАЮЩИЙ РАЧЕТ
// ============================================================================

class RatchetFlashing : public StochasticTask {
public:
    RatchetFlashing(BaseDiffusionSolver& solver,
                   Plotter* plotter,
                   unsigned short task_id,
                   const std::string& task_name)
        : StochasticTask(solver, plotter, task_id, task_name) {}

protected:
    std::unique_ptr<Potential1D> create_potential(const DiffusionParameters& params) override {
        return std::make_unique<RatchetPotential>(params.ratchet_amplitude,
                                                  params.ratchet_wavenumber);
    }

    void visualize(const DiffusionResult& result,
                  const DiffusionParameters& params) override {
        RatchetPotential pot(params.ratchet_amplitude, params.ratchet_wavenumber);
        visualizer_.plot_potential(pot, params, "Asymmetric ratchet potential");
        visualizer_.plot_trajectories(result, 10, "Trajectories (flashing ratchet)");
        visualizer_.plot_mean_position(result, "Mean position (RATCHET EFFECT!)");
        visualizer_.plot_statistics(result, "Flashing Ratchet");
    }
};

// ============================================================================
// ЗАДАЧА 5: НАКЛОННЫЙ РАЧЕТ
// ============================================================================

class RatchetTilting : public StochasticTask {
public:
    RatchetTilting(BaseDiffusionSolver& solver,
                  Plotter* plotter,
                  unsigned short task_id,
                  const std::string& task_name)
        : StochasticTask(solver, plotter, task_id, task_name) {}

protected:
    std::unique_ptr<Potential1D> create_potential(const DiffusionParameters& params) override {
        return std::make_unique<RatchetPotential>(params.ratchet_amplitude,
                                                  params.ratchet_wavenumber);
    }

    void visualize(const DiffusionResult& result,
                  const DiffusionParameters& params) override {
        visualizer_.plot_trajectories(result, 10, "Trajectories (tilting ratchet)");
        visualizer_.plot_mean_position(result, "Mean position (tilting)");
        visualizer_.plot_statistics(result, "Tilting Ratchet");
    }
};

} // namespace special

#endif // NUMERICAL_METHODS_IN_PHYSICS_STOCHASTIC_TASKS_H