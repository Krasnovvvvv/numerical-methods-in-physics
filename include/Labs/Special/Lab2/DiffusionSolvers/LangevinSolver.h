#ifndef NUMERICAL_METHODS_IN_PHYSICS_LANGEVIN_SOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_LANGEVIN_SOLVER_H

#pragma once

#include "Labs/Special/Lab2/Base/BaseDiffusionSolver.h"
#include "Labs/Special/Lab2/Base/DiffusionParameters.h"
#include "Labs/Special/Lab2/Base/Potential1D.h"

#include <random>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace special {

/**
 * @class LangevinSolver
 * @brief Решатель уравнения Ланжевена (1D) в режиме сверхдемпфирования и инерциальном режиме
 *
 * Реализует схему Эйлера–Маруямы для стохастического дифференциального уравнения:
 *
 * **Сверхдемпфирвованный режим (overdamped = true):**
 * \f[
 *   dx = \frac{D}{k_B T} F(x, t) \, dt + \sqrt{2D} \, dW
 * \f]
 *
 * **Инерциальный режим (overdamped = false):**
 * \f[
 *   dx = v \, dt
 * \f]
 * \f[
 *   dv = \left( \frac{F(x, t)}{m} - \frac{\gamma v}{m} \right) dt + \sqrt{\frac{2\gamma k_B T}{m^2}} \, dW
 * \f]
 *
 * Вычисляет все 4 момента: ⟨x⟩, ⟨x²⟩, ⟨x³⟩, ⟨x⁴⟩ для анализа качества метода.
 */
class LangevinSolver : public BaseDiffusionSolver {
public:
    /**
     * @brief Конструктор
     * @param seed Seed для генератора случайных чисел (для воспроизводимости)
     */
    explicit LangevinSolver(unsigned int seed = 42)
        : rng_(seed),
          normal_dist_(0.0, 1.0) {}

    ~LangevinSolver() override = default;

    /**
     * @return Название метода решения
     */
    std::string name() const override {
        return "LangevinSolver (Euler-Maruyama)";
    }

    /**
     * @brief Основной метод решения
     *
     * Численно интегрирует уравнение Ланжевена. Вычисляет траектории всех частиц
     * и все 4 момента распределения на каждом временном шаге.
     *
     * @param params Параметры моделирования
     * @param potential Потенциальная функция U(x, t)
     * @return DiffusionResult с траекториями, моментами и статистикой
     * @throws std::invalid_argument если параметры невалидны
     */
    DiffusionResult solve(const DiffusionParameters& params,
                          const Potential1D& potential) override {
        if (params.n_particles <= 0 || params.n_steps <= 0) {
            throw std::invalid_argument("n_particles and n_steps must be > 0");
        }
        if (params.dt <= 0.0) {
            throw std::invalid_argument("dt must be > 0");
        }

        DiffusionResult result;
        const int n_particles = params.n_particles;
        const int n_steps = params.n_steps;
        const double dt = params.dt;

        // Инициализация результатов
        result.t.resize(n_steps + 1);
        result.mean_x.resize(n_steps + 1, 0.0);
        result.mean_x2.resize(n_steps + 1, 0.0);
        result.mean_x3.resize(n_steps + 1, 0.0);  // НОВОЕ
        result.mean_x4.resize(n_steps + 1, 0.0);  // НОВОЕ
        result.displacement_squared.resize(n_steps + 1, 0.0);
        result.trajectories.assign(n_particles, std::vector<double>(n_steps + 1, 0.0));

        // Начальные условия
        std::vector<double> x(n_particles, params.x0);
        std::vector<double> v(n_particles, 0.0);

        // Интегрирование
        for (int i = 0; i <= n_steps; ++i) {
            double t = i * dt;
            result.t[i] = t;

            double sum_x = 0.0;
            double sum_x2 = 0.0;
            double sum_x3 = 0.0;
            double sum_x4 = 0.0;
            double sum_disp2 = 0.0;

            // Накопление статистики
            for (int p = 0; p < n_particles; ++p) {
                double dx = x[p] - params.x0;
                double x_p = x[p];

                sum_x += x_p;
                sum_x2 += x_p * x_p;
                sum_x3 += x_p * x_p * x_p;
                sum_x4 += x_p * x_p * x_p * x_p;
                sum_disp2 += dx * dx;

                result.trajectories[p][i] = x_p;
            }

            double invN = 1.0 / static_cast<double>(n_particles);
            result.mean_x[i] = sum_x * invN;
            result.mean_x2[i] = sum_x2 * invN;
            result.mean_x3[i] = sum_x3 * invN;  // НОВОЕ
            result.mean_x4[i] = sum_x4 * invN;  // НОВОЕ
            result.displacement_squared[i] = sum_disp2 * invN;

            if (i == n_steps) {
                break;
            }

            // Интегрирование на один шаг
            if (params.overdamped) {
                step_overdamped(x, t, dt, params, potential);
            } else {
                step_inertial(x, v, t, dt, params, potential);
            }

            // Применение граничных условий
            if (params.use_periodic_bc) {
                apply_periodic_bc(x, params);
            }
        }

        result.steps = n_steps;
        return result;
    }

    /**
     * @brief Один шаг для анимации (переопределение виртуального метода)
     */
    std::vector<double> step(const DiffusionParameters& params,
                             const Potential1D& potential,
                             const std::vector<double>& state,
                             double time) override {
        std::vector<double> x = state;
        std::vector<double> v(state.size(), 0.0);

        if (params.overdamped) {
            step_overdamped(x, time, params.dt, params, potential);
        } else {
            step_inertial(x, v, time, params.dt, params, potential);
        }

        if (params.use_periodic_bc) {
            apply_periodic_bc(x, params);
        }

        return x;
    }

private:
    std::mt19937 rng_;
    std::normal_distribution<double> normal_dist_;

    /**
     * @brief Один шаг в сверхдемпфированном режиме
     */
    void step_overdamped(std::vector<double>& x,
                         double t,
                         double dt,
                         const DiffusionParameters& params,
                         const Potential1D& potential) {
        const int n_particles = static_cast<int>(x.size());
        const double D = params.diffusion_coeff;
        const double kBT = params.kB_T;
        const double sqrt_2Ddt = std::sqrt(2.0 * D * dt);

        for (int p = 0; p < n_particles; ++p) {
            double F = potential.F(x[p], t);
            double drift = (D / kBT) * F * dt;
            double noise = sqrt_2Ddt * normal_dist_(rng_);
            x[p] += drift + noise;
        }
    }

    /**
     * @brief Один шаг в инерциальном режиме
     */
    void step_inertial(std::vector<double>& x,
                       std::vector<double>& v,
                       double t,
                       double dt,
                       const DiffusionParameters& params,
                       const Potential1D& potential) {
        const int n_particles = static_cast<int>(x.size());
        const double m = params.mass;
        const double gamma = params.friction;
        const double kBT = params.kB_T;

        const double sqrt_2gamma_kBT_dt_over_m2 =
            std::sqrt(2.0 * gamma * kBT * dt / (m * m));

        for (int p = 0; p < n_particles; ++p) {
            double F = potential.F(x[p], t);
            double noise = sqrt_2gamma_kBT_dt_over_m2 * normal_dist_(rng_);

            v[p] += (F / m - gamma * v[p] / m) * dt + noise;
            x[p] += v[p] * dt;
        }
    }

    /**
     * @brief Применение периодических граничных условий
     */
    void apply_periodic_bc(std::vector<double>& x,
                           const DiffusionParameters& params) {
        if (!params.use_periodic_bc) return;

        const double x_min = params.x_min;
        const double x_max = params.x_max;
        const double len = x_max - x_min;

        if (len <= 0.0) return;

        for (double& xi : x) {
            while (xi < x_min) xi += len;
            while (xi > x_max) xi -= len;
        }
    }
};

} // namespace special

#endif // NUMERICAL_METHODS_IN_PHYSICS_LANGEVIN_SOLVER_H