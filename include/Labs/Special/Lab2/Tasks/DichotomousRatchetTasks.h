#ifndef DICHOTOMIC_RATCHET_TASKS_H
#define DICHOTOMIC_RATCHET_TASKS_H

#pragma once

#include "Labs/Special/Lab2/Base/StochasticTask.h"
#include "Labs/Special/Lab2/Base/DichotomousRatchetPotential.h"
#include "Labs/Special/Lab2/Base/DichotomousNoiseGenerator.h"
#include "Labs/Special/Lab2/DiffusionSolvers/LangevinSolver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>

namespace special {

// ============================================================================
// ЧАСТЬ 1: ЗАДАНИЯ ПО ГЕНЕРИРОВАНИЮ ДИХОТОМНОГО ШУМА
// ============================================================================

/**
 * @brief ЗАДАНИЕ 1: Генерирование дихотомного шума и визуализация
 */
class NoiseGenerationTask : public StochasticTask {
public:
    NoiseGenerationTask(BaseDiffusionSolver& solver,
                       Plotter* plotter,
                       unsigned short task_id,
                       const std::string& task_name)
        : StochasticTask(solver, plotter, task_id, task_name)
    {
    }

    void run(const DiffusionParameters& params) {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "PART 1, TASK 1: Dichotomic Noise Generation\n";
        std::cout << std::string(70, '=') << "\n\n";

        // ЗАДАНИЕ 1: симметричный шум
        std::cout << "PART 1, SUBTASK 1a): Symmetric noise (γ_a = γ_b = Γ/2)\n";
        std::cout << "Parameters: a = 1, b = -1, τ_C = 0.05, 0.4\n\n";

        run_symmetric_noise(0.05);
        run_symmetric_noise(0.4);

        // ЗАДАНИЕ 2: асимметричный шум
        std::cout << "\n\nPART 1, SUBTASK 2): Asymmetric noise (a = -b, γ_a ≠ γ_b)\n";
        std::cout << "Parameters: a = 2, b = -1, τ_C = 0.05, 0.4\n\n";

        run_asymmetric_noise(0.05);
        run_asymmetric_noise(0.4);

        std::cout << "\n✓ All noise generation tasks completed\n";
    }

protected:
    void run_symmetric_noise(double tau_c) {
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        std::cout << "τ_C = " << std::fixed << std::setprecision(2) << tau_c << "\n";

        double Gamma = 1.0 / tau_c;
        double gamma_a = Gamma / 2.0;
        double gamma_b = Gamma / 2.0;

        DichotomousNoiseGenerator gen(1.0, -1.0, gamma_a, gamma_b);

        double dt = tau_c / 20.0;
        int n_steps = static_cast<int>(10.0 / dt);

        std::cout << "Generated: dt = " << dt << ", n_steps = " << n_steps;

        auto sigma = gen.generate(dt, n_steps);

        double mean = gen.mean(sigma);
        double m2 = gen.second_moment(sigma);

        std::cout << "\nMoments:\n";
        std::cout << "  ⟨σ⟩ = " << std::scientific << mean << " (should be ≈ 0)\n";
        std::cout << "  ⟨σ²⟩ = " << m2 << " (should be ≈ 1.0)\n";

        bool mean_ok = std::abs(mean) < 1e-3;
        bool m2_ok = std::abs(m2 - 1.0) < 1e-2;

        if (mean_ok && m2_ok) {
            std::cout << "✓ PASS: Noise properties are correct\n";
        } else {
            std::cout << "✗ FAIL: Check noise generation!\n";
        }

        std::string filename = "noise_symmetric_tau" + std::to_string(static_cast<int>(tau_c*100)) + ".txt";
        std::ofstream fout(filename);
        for (int i = 0; i < std::min(500, static_cast<int>(sigma.size())); ++i) {
            fout << i * dt << " " << sigma[i] << "\n";
        }
        fout.close();
        std::cout << "✓ Saved to: " << filename << "\n";

        std::cout << "\n";
    }

    void run_asymmetric_noise(double tau_c) {
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        std::cout << "τ_C = " << std::fixed << std::setprecision(2) << tau_c << "\n";

        double Gamma = 1.0 / tau_c;
        double gamma_b = Gamma / 3.0;
        double gamma_a = 2.0 * Gamma / 3.0;

        DichotomousNoiseGenerator gen(2.0, -1.0, gamma_a, gamma_b);

        double dt = tau_c / 20.0;
        int n_steps = static_cast<int>(10.0 / dt);

        std::cout << "Parameters: a = 2.0, b = -1.0\n";
        std::cout << "  γ_a = " << gamma_a << ", γ_b = " << gamma_b << "\n";
        std::cout << "Generated: dt = " << dt << ", n_steps = " << n_steps;

        auto sigma = gen.generate(dt, n_steps);

        double mean = gen.mean(sigma);
        double m2 = gen.second_moment(sigma);

        std::cout << "\nMoments:\n";
        std::cout << "  ⟨σ⟩ = " << std::scientific << mean << " (should be ≈ 0)\n";
        std::cout << "  ⟨σ²⟩ = " << m2 << "\n";

        bool mean_ok = std::abs(mean) < 1e-3;

        if (mean_ok) {
            std::cout << "✓ PASS: Zero mean condition satisfied\n";
        } else {
            std::cout << "✗ FAIL: Mean is not zero!\n";
        }

        std::string filename = "noise_asymmetric_tau" + std::to_string(static_cast<int>(tau_c*100)) + ".txt";
        std::ofstream fout(filename);
        for (int i = 0; i < std::min(500, static_cast<int>(sigma.size())); ++i) {
            fout << i * dt << " " << sigma[i] << "\n";
        }
        fout.close();
        std::cout << "✓ Saved to: " << filename << "\n";

        std::cout << "\n";
    }

    std::unique_ptr<Potential1D> create_potential(
        const DiffusionParameters& params) override
    {
        return nullptr;
    }
};

// ============================================================================
// ЗАДАНИЕ 3: Проверка алгоритма (вычисление моментов)
// ============================================================================

class NoiseVerificationTask : public StochasticTask {
public:
    NoiseVerificationTask(BaseDiffusionSolver& solver,
                         Plotter* plotter,
                         unsigned short task_id,
                         const std::string& task_name)
        : StochasticTask(solver, plotter, task_id, task_name)
    {
    }

    void run(const DiffusionParameters& params) {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "PART 1, TASK 3: Noise Verification (First Moment)\n";
        std::cout << std::string(70, '=') << "\n\n";

        std::cout << "Check that ⟨σ(t)⟩ is close to zero for different configurations\n\n";

        // Тест 1: симметричный с a=1, b=-1
        {
            std::cout << "Test 1: a=1, b=-1, γ_a=1.0, γ_b=1.0\n";
            DichotomousNoiseGenerator gen(1.0, -1.0, 1.0, 1.0);
            auto sigma = gen.generate(0.01, 10000);
            double mean = gen.mean(sigma);
            std::cout << "  ⟨σ⟩ = " << std::scientific << mean;
            if (std::abs(mean) < 1e-3) {
                std::cout << " ✓ PASS\n";
            } else {
                std::cout << " ✗ FAIL\n";
            }
        }

        // Тест 2: асимметричный с a=2, b=-1
        {
            std::cout << "\nTest 2: a=2, b=-1, γ_a=2/3, γ_b=1/3 (for zero mean)\n";
            DichotomousNoiseGenerator gen(2.0, -1.0, 2.0/3.0, 1.0/3.0);
            auto sigma = gen.generate(0.01, 10000);
            double mean = gen.mean(sigma);
            std::cout << "  ⟨σ⟩ = " << std::scientific << mean;
            if (std::abs(mean) < 1e-3) {
                std::cout << " ✓ PASS\n";
            } else {
                std::cout << " ✗ FAIL\n";
            }
        }

        std::cout << "\n✓ Verification complete\n";
    }

protected:
    std::unique_ptr<Potential1D> create_potential(
        const DiffusionParameters& params) override
    {
        return nullptr;
    }
};

// ============================================================================
// ЧАСТЬ 2: ЗАДАНИЯ ПО РАЧЕТ-ПОТЕНЦИАЛАМ
// ============================================================================

/**
 * @class DichotomicLangevinTaskBase
 * @brief Базовый класс для задач Ланжевена с дихотомным шумом
 */
class DichotomicLangevinTaskBase : public StochasticTask {
protected:
    DichotomicLangevinTaskBase(BaseDiffusionSolver& solver,
                               Plotter* plotter,
                               unsigned short task_id,
                               const std::string& task_name)
        : StochasticTask(solver, plotter, task_id, task_name)
    {
    }

    /**
     * @brief Решить уравнение Ланжевена с дихотомным шумом
     */
    DiffusionResult solve_with_dichotomic_noise(
        const DiffusionParameters& params,
        DichotomousRatchetPotential* ratchet_potential)
    {
        if (!ratchet_potential) {
            throw std::invalid_argument("ratchet_potential is null");
        }

        auto generator = ratchet_potential->generator();
        int n_noise_steps = params.n_steps + 1;
        auto sigma_values = generator->generate(params.dt, n_noise_steps);

        DiffusionResult result;
        result.t.resize(params.n_steps + 1);
        result.mean_x.resize(params.n_steps + 1, 0.0);
        result.mean_x2.resize(params.n_steps + 1, 0.0);
        result.mean_x3.resize(params.n_steps + 1, 0.0);
        result.mean_x4.resize(params.n_steps + 1, 0.0);
        result.displacement_squared.resize(params.n_steps + 1, 0.0);
        result.trajectories.assign(params.n_particles,
                                  std::vector<double>(params.n_steps + 1, 0.0));

        std::vector<double> x(params.n_particles, params.x0);

        for (int i = 0; i <= params.n_steps; ++i) {
            double t = i * params.dt;
            result.t[i] = t;

            if (i < static_cast<int>(sigma_values.size())) {
                ratchet_potential->set_sigma(sigma_values[i], t);
            }

            double sum_x = 0.0, sum_x2 = 0.0, sum_x3 = 0.0, sum_x4 = 0.0;
            double sum_disp2 = 0.0;

            for (int p = 0; p < params.n_particles; ++p) {
                double xp = x[p];
                double dx = xp - params.x0;

                sum_x += xp;
                sum_x2 += xp * xp;
                sum_x3 += xp * xp * xp;
                sum_x4 += xp * xp * xp * xp;
                sum_disp2 += dx * dx;

                result.trajectories[p][i] = xp;
            }

            double invN = 1.0 / params.n_particles;
            result.mean_x[i] = sum_x * invN;
            result.mean_x2[i] = sum_x2 * invN;
            result.mean_x3[i] = sum_x3 * invN;
            result.mean_x4[i] = sum_x4 * invN;
            result.displacement_squared[i] = sum_disp2 * invN;

            if (i < params.n_steps) {
                x = solver_.step(params, *ratchet_potential, x, t);

                if (params.use_periodic_bc) {
                    const double len = params.x_max - params.x_min;
                    for (double& xi : x) {
                        while (xi < params.x_min) xi += len;
                        while (xi > params.x_max) xi -= len;
                    }
                }
            }
        }

        result.steps = params.n_steps;
        return result;
    }

    /**
     * @brief Вычислить среднюю скорость рачета
     */
    double compute_ratchet_velocity(const DiffusionResult& result,
                                    const DiffusionParameters& params) const
    {
        if (result.mean_x.empty() || result.t.empty()) {
            return 0.0;
        }

        double x_final = result.mean_x.back();
        double x_initial = result.mean_x.front();
        double t_max = result.t.back();

        if (t_max > 1e-10) {
            return (x_final - x_initial) / t_max;
        }
        return 0.0;
    }
};

// ============================================================================
// ЗАДАНИЕ 1A: On-Off рачет (симметричный шум)
// ============================================================================

class OnOffRatchetSymmetric : public DichotomicLangevinTaskBase {
public:
    OnOffRatchetSymmetric(BaseDiffusionSolver& solver,
                         Plotter* plotter,
                         unsigned short task_id,
                         const std::string& task_name,
                         double V0, double L, double tau_c)
        : DichotomicLangevinTaskBase(solver, plotter, task_id, task_name),
          V0_(V0), L_(L), tau_c_(tau_c)
    {
    }

    void run(const DiffusionParameters& params) {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "PART 2, TASK 1A: On-Off Ratchet (Symmetric Noise)\n";
        std::cout << "V₀=" << V0_ << ", L=" << L_ << ", τ_C=" << tau_c_ << "\n";
        std::cout << std::string(70, '=') << "\n\n";

        auto potential = create_potential(params);
        auto ratchet = dynamic_cast<DichotomousRatchetPotential*>(potential.get());
        if (!ratchet) throw std::runtime_error("Failed to create potential");

        DiffusionResult result = solve_with_dichotomic_noise(params, ratchet);

        double v = compute_ratchet_velocity(result, params);
        std::cout << "Mean velocity: v = " << std::scientific << v << "\n";
        std::cout << "Expected: v ≈ 0 (symmetric noise → no drift)\n";
        std::cout << std::defaultfloat;

        visualize(result, params);
        postprocess(result, params);
    }

protected:
    std::unique_ptr<Potential1D> create_potential(
        const DiffusionParameters& params) override
    {
        double Gamma = 1.0 / tau_c_;
        double gamma_a = Gamma / 2.0;
        double gamma_b = Gamma / 2.0;

        return std::make_unique<DichotomousRatchetPotential>(
            V0_, L_, 0.25,
            gamma_a, gamma_b,
            DichotomousRatchetPotential::ModulationType::ON_OFF,
            0.0);
    }

    void visualize(const DiffusionResult& result,
                  const DiffusionParameters& params)
    {
        std::cout << "\n✓ Generating plots...\n";
        visualizer_.plot_trajectories(result, 5,
            "Task 1A: On-Off Ratchet - Sample Trajectories");
        visualizer_.plot_mean_position(result,
            "Task 1A: On-Off Ratchet - Mean Position ⟨x(t)⟩");
        visualizer_.plot_statistics(result,
            "Task 1A: On-Off Ratchet - Moments");
    }

    void postprocess(const DiffusionResult& result,
                    const DiffusionParameters& params)
    {
        std::cout << "Statistics:\n";
        std::cout << "  Final ⟨x⟩ = " << result.mean_x.back() << "\n";
        std::cout << "  Displacement² = " << result.displacement_squared.back() << "\n";
    }

private:
    double V0_, L_, tau_c_;
};

// ============================================================================
// ЗАДАНИЕ 2B: Flipping рачет (асимметричный шум)
// ============================================================================

class FlippingRatchetAsymmetric : public DichotomicLangevinTaskBase {
public:
    FlippingRatchetAsymmetric(BaseDiffusionSolver& solver,
                             Plotter* plotter,
                             unsigned short task_id,
                             const std::string& task_name,
                             double V0, double L, double tau_c,
                             double ratio_gamma)
        : DichotomicLangevinTaskBase(solver, plotter, task_id, task_name),
          V0_(V0), L_(L), tau_c_(tau_c), ratio_gamma_(ratio_gamma)
    {
    }

    void run(const DiffusionParameters& params) {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "PART 2, TASK 2B: Flipping Ratchet (Asymmetric Noise)\n";
        std::cout << "V₀=" << V0_ << ", L=" << L_ << ", τ_C=" << tau_c_
                 << ", γ_a/γ_b=" << ratio_gamma_ << "\n";
        std::cout << std::string(70, '=') << "\n\n";

        auto potential = create_potential(params);
        auto ratchet = dynamic_cast<DichotomousRatchetPotential*>(potential.get());
        if (!ratchet) throw std::runtime_error("Failed to create potential");

        DiffusionResult result = solve_with_dichotomic_noise(params, ratchet);

        double v = compute_ratchet_velocity(result, params);
        std::cout << "Mean velocity: v = " << std::scientific << v << "\n";
        std::cout << "Expected: v ≠ 0 (asymmetric noise → ratchet effect)\n";
        std::cout << std::defaultfloat;

        visualize(result, params);
        postprocess(result, params);
    }

protected:
    std::unique_ptr<Potential1D> create_potential(
        const DiffusionParameters& params) override
    {
        double Gamma = 1.0 / tau_c_;
        double gamma_a = Gamma * ratio_gamma_ / (1.0 + ratio_gamma_);
        double gamma_b = Gamma / (1.0 + ratio_gamma_);

        return std::make_unique<DichotomousRatchetPotential>(
            V0_, L_, 0.25,
            gamma_a, gamma_b,
            DichotomousRatchetPotential::ModulationType::FLIPPING,
            0.0);
    }

    void visualize(const DiffusionResult& result,
                  const DiffusionParameters& params)
    {
        std::cout << "\n✓ Generating plots...\n";
        visualizer_.plot_trajectories(result, 5,
            "Task 2B: Flipping Ratchet - Sample Trajectories");
        visualizer_.plot_mean_position(result,
            "Task 2B: Flipping Ratchet - Mean Position ⟨x(t)⟩");
        visualizer_.plot_statistics(result,
            "Task 2B: Flipping Ratchet - Moments");
    }

    void postprocess(const DiffusionResult& result,
                    const DiffusionParameters& params)
    {
        std::cout << "Statistics:\n";
        std::cout << "  Final ⟨x⟩ = " << result.mean_x.back() << "\n";
        std::cout << "  Displacement² = " << result.displacement_squared.back() << "\n";
    }

private:
    double V0_, L_, tau_c_;
    double ratio_gamma_;
};

// ============================================================================
// ЗАДАНИЕ 2V*: Зависимость скорости от τ_C
// ============================================================================

class RatchetVelocityVsTauC : public DichotomicLangevinTaskBase {
public:
    RatchetVelocityVsTauC(BaseDiffusionSolver& solver,
                         Plotter* plotter,
                         unsigned short task_id,
                         const std::string& task_name,
                         double V0, double L, double ratio_gamma)
        : DichotomicLangevinTaskBase(solver, plotter, task_id, task_name),
          V0_(V0), L_(L), ratio_gamma_(ratio_gamma)
    {
    }

    void run(const DiffusionParameters& params_base) {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "PART 2, TASK 2V*: Ratchet Velocity vs τ_C\n";
        std::cout << "Scanning τ_C for optimal velocity (Bell-shaped curve)\n";
        std::cout << std::string(70, '=') << "\n\n";

        std::vector<double> tau_c_values;
        std::vector<double> velocities;

        double tau_c_min = 0.01;
        double tau_c_max = 5.0;
        int n_points = 20;

        for (int i = 0; i < n_points; ++i) {
            double tau_c = tau_c_min + (tau_c_max - tau_c_min) * i / (n_points - 1);
            tau_c_values.push_back(tau_c);

            DiffusionParameters params = params_base;
            params.dt = tau_c / 20.0;
            params.n_steps = static_cast<int>(10.0 / params.dt);

            std::cout << "τ_C = " << std::fixed << std::setprecision(3) << tau_c
                     << ", dt = " << params.dt << ", steps = " << params.n_steps
                     << " ... ";
            std::cout.flush();

            double Gamma = 1.0 / tau_c;
            double gamma_a = Gamma * ratio_gamma_ / (1.0 + ratio_gamma_);
            double gamma_b = Gamma / (1.0 + ratio_gamma_);

            auto potential = std::make_unique<DichotomousRatchetPotential>(
                V0_, L_, 0.25, gamma_a, gamma_b,
                DichotomousRatchetPotential::ModulationType::FLIPPING, 0.0);

            try {
                auto ratchet = dynamic_cast<DichotomousRatchetPotential*>(potential.get());
                DiffusionResult result = solve_with_dichotomic_noise(params, ratchet);
                double v = compute_ratchet_velocity(result, params);
                velocities.push_back(v);
                std::cout << "v = " << std::scientific << v << "\n";
            } catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << "\n";
                velocities.push_back(0.0);
            }
        }

        std::ofstream fout("ratchet_velocity_vs_tau_c.txt");
        for (size_t i = 0; i < tau_c_values.size(); ++i) {
            fout << std::scientific << tau_c_values[i] << " " << velocities[i] << "\n";
        }
        fout.close();

        std::cout << "\n✓ Results saved to: ratchet_velocity_vs_tau_c.txt\n";
        std::cout << "Expected: Bell-shaped curve with maximum\n";
    }

protected:
    std::unique_ptr<Potential1D> create_potential(
        const DiffusionParameters& params) override
    {
        double Gamma = 1.0;
        double gamma_a = Gamma * ratio_gamma_ / (1.0 + ratio_gamma_);
        double gamma_b = Gamma / (1.0 + ratio_gamma_);

        return std::make_unique<DichotomousRatchetPotential>(
            V0_, L_, 0.25, gamma_a, gamma_b,
            DichotomousRatchetPotential::ModulationType::FLIPPING, 0.0);
    }

private:
    double V0_, L_;
    double ratio_gamma_;
};

} // namespace special

#endif // DICHOTOMIC_RATCHET_TASKS_H