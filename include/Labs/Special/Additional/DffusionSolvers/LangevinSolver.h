#ifndef NUMERICAL_METHODS_IN_PHYSICS_LANGEVINSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_LANGEVINSOLVER_H

#pragma once

#include <vector>
#include <random>
#include <cmath>
#include <functional>
#include <stdexcept>
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"

// ============================================================================
// Структура для одной траектории (старая версия, для совместимости)
// ============================================================================
struct LangevinTrajectory {
    std::vector<double> t;           // время
    std::vector<double> x;           // позиция
    double L;                        // период
    double mean_position;
    double mean_velocity;
    double diffusion_coeff;
};

// ============================================================================
// Структура для ансамбля частиц
// ============================================================================
struct LangevinEnsembleResult {
    std::vector<double> t;           // временные точки
    std::vector<double> mean_x;      // средняя позиция по ансамблю <x(t)>
    std::vector<double> mean_x2;     // <x²(t)> для анализа
    double L;                        // период
    double mean_velocity;            // средняя скорость по ансамблю
    double diffusion_coeff;          // коэффициент диффузии
    std::size_t n_particles;         // число частиц
};

// ============================================================================
// Решатель Ланжевена
// ============================================================================
class LangevinSolver {
public:
    using DVDXFunc = std::function<double(double)>;

    enum class ModulationType {
        CONSTANT,              // f(t) = 1.0
        DICHOTOM_SYMMETRIC,    // f(t) = (1/2)(1 + σ(t))
        EPSILON_PLUS_DICHOTOM  // f(t) = ε + σ(t)
    };

    LangevinSolver(DVDXFunc dVdx,
                   double V0,
                   double L,
                   double dt,
                   ModulationType mod_type = ModulationType::CONSTANT,
                   unsigned int seed = std::random_device{}())
        : dVdx_(dVdx),
          V0_(V0),
          L_(L),
          dt_(dt),
          mod_type_(mod_type),
          rng_(seed),
          gaussian_(0.0, 1.0)
    {
        if (dt_ <= 0.0) throw std::invalid_argument("dt must be > 0");
        if (L_ <= 0.0) throw std::invalid_argument("L must be > 0");
        if (V0_ < 0.0) throw std::invalid_argument("V0 must be >= 0");
        sqrt2_coeff_ = std::sqrt(2.0);
    }

    // ========================================================================
    // Старый метод: одна частица
    // ========================================================================
    LangevinTrajectory solve(DichotomicNoise& noise,
                             std::size_t N,
                             double x0 = 0.0,
                             double epsilon = 0.0,
                             std::size_t burn_in = 0)
    {
        auto dichotom_profile = noise.generate(N, 0, true);

        LangevinTrajectory res;
        res.t.resize(N);
        res.x.resize(N);
        res.L = L_;

        double x = std::fmod(x0, L_);
        if (x < 0.0) x += L_;

        for (std::size_t n = 0; n < N; ++n) {
            res.t[n] = static_cast<double>(n) * dt_;
            res.x[n] = x;

            double sigma_n = dichotom_profile.s[n];
            double f_t = compute_f_t(sigma_n, epsilon);

            double zeta_n = gaussian_(rng_);
            double dVdx_n = dVdx_(x);
            double dxdt = -V0_ * dVdx_n * f_t + sqrt2_coeff_ * zeta_n;
            x += dxdt * dt_;

            x = std::fmod(x, L_);
            if (x < 0.0) x += L_;
        }

        if (burn_in > N) burn_in = N;

        // средняя позиция
        double sum_x = 0.0;
        std::size_t cnt = 0;
        for (std::size_t n = burn_in; n < N; ++n) {
            sum_x += res.x[n];
            ++cnt;
        }
        res.mean_position = (cnt > 0) ? (sum_x / static_cast<double>(cnt)) : 0.0;

        // средняя скорость
        double t_max = (N - burn_in) * dt_;
        if (N > burn_in + 1 && t_max > 0.0) {
            double sum_disp = 0.0;
            double x0_eff = res.x[burn_in];
            for (std::size_t n = burn_in; n < N; ++n) {
                double dx = res.x[n] - x0_eff;
                if (dx > L_ / 2.0) dx -= L_;
                if (dx < -L_ / 2.0) dx += L_;
                sum_disp += dx;
            }
            std::size_t N_eff = N - burn_in;
            res.mean_velocity = (1.0 / t_max) * (1.0 / static_cast<double>(N_eff)) * sum_disp;
        } else {
            res.mean_velocity = 0.0;
        }

        // коэффициент диффузии
        if (N > burn_in + 1) {
            double sum_dx2 = 0.0;
            std::size_t cnt_dx = 0;
            for (std::size_t n = burn_in + 1; n < N; ++n) {
                double dx = res.x[n] - res.x[n - 1];
                if (dx > L_ / 2.0) dx -= L_;
                if (dx < -L_ / 2.0) dx += L_;
                sum_dx2 += dx * dx;
                ++cnt_dx;
            }
            double mean_dx2 = sum_dx2 / static_cast<double>(cnt_dx);
            res.diffusion_coeff = mean_dx2 / (2.0 * dt_);
        } else {
            res.diffusion_coeff = 0.0;
        }

        return res;
    }

    // ========================================================================
    // НОВЫЙ метод: ансамбль из n_particles частиц в ОДНОМ общем шуме
    // ========================================================================
    LangevinEnsembleResult solve_ensemble(DichotomicNoise& noise,
                                          std::size_t N,
                                          std::size_t n_particles,
                                          double x0 = 0.0,
                                          double epsilon = 0.0,
                                          std::size_t burn_in = 0)
    {

        auto dichotom_profile = noise.generate(N, 0, true);

        LangevinEnsembleResult res;
        res.t.resize(N);
        res.mean_x.resize(N);
        res.mean_x2.resize(N);
        res.L = L_;
        res.n_particles = n_particles;

        // Массив частиц: x[p] — позиция частицы p
        std::vector<double> x(n_particles, 0.0);
        for (auto& xi : x) {
            xi = std::fmod(x0, L_);
            if (xi < 0.0) xi += L_;
        }

        // Интегрирование: все частицы в общем шуме σ(t)
        for (std::size_t n = 0; n < N; ++n) {
            res.t[n] = static_cast<double>(n) * dt_;

            double sigma_n = dichotom_profile.s[n];
            double f_t = compute_f_t(sigma_n, epsilon);

            // Сумма координат и квадратов координат
            double sum_x = 0.0;
            double sum_x2 = 0.0;

            // Шаг интегрирования для каждой частицы
            for (std::size_t p = 0; p < n_particles; ++p) {
                // Каждой частице свой гауссовский шум
                double zeta_n = gaussian_(rng_);
                double dVdx_n = dVdx_(x[p]);
                double dxdt = -V0_ * dVdx_n * f_t + sqrt2_coeff_ * zeta_n;
                x[p] += dxdt * dt_;

                // Периодические ГУ
                x[p] = std::fmod(x[p], L_);
                if (x[p] < 0.0) x[p] += L_;

                sum_x += x[p];
                sum_x2 += x[p] * x[p];
            }

            res.mean_x[n] = sum_x / static_cast<double>(n_particles);
            res.mean_x2[n] = sum_x2 / static_cast<double>(n_particles);
        }

        if (burn_in > N) burn_in = N;

        // === Средняя скорость по ансамблю (после burn_in) ===
        double t_max = (N - burn_in) * dt_;
        if (N > burn_in + 1 && t_max > 0.0) {
            double sum_all_disp = 0.0;
            double x0_eff_mean = res.mean_x[burn_in];

            for (std::size_t n = burn_in; n < N; ++n) {
                double dx = res.mean_x[n] - x0_eff_mean;
                if (dx > L_ / 2.0) dx -= L_;
                if (dx < -L_ / 2.0) dx += L_;
                sum_all_disp += dx;
            }

            std::size_t N_eff = N - burn_in;
            res.mean_velocity = (1.0 / t_max) * (1.0 / static_cast<double>(N_eff)) * sum_all_disp;
        } else {
            res.mean_velocity = 0.0;
        }

        // === Коэффициент диффузии (по MSD после burn_in) ===
        if (N > burn_in + 1) {
            double sum_msd = 0.0;
            std::size_t cnt = 0;

            for (std::size_t n = burn_in + 1; n < N; ++n) {
                double dmean_x = res.mean_x[n] - res.mean_x[n - 1];
                if (dmean_x > L_ / 2.0) dmean_x -= L_;
                if (dmean_x < -L_ / 2.0) dmean_x += L_;

                sum_msd += dmean_x * dmean_x;
                ++cnt;
            }

            double mean_msd = sum_msd / static_cast<double>(cnt);
            res.diffusion_coeff = mean_msd / (2.0 * dt_);
        } else {
            res.diffusion_coeff = 0.0;
        }

        return res;
    }

    // геттеры
    double L() const { return L_; }
    double dt() const { return dt_; }
    double V0() const { return V0_; }

private:
    DVDXFunc dVdx_;
    double V0_, L_, dt_;
    ModulationType mod_type_;
    double sqrt2_coeff_;
    std::mt19937 rng_;
    std::normal_distribution<double> gaussian_;

    // Вспомогательная функция для вычисления модуляции f(t)
    double compute_f_t(double sigma_n, double epsilon)
    {
        switch (mod_type_) {
            case ModulationType::CONSTANT:
                return 1.0;
            case ModulationType::DICHOTOM_SYMMETRIC:
                return 0.5 * (1.0 + sigma_n);
            case ModulationType::EPSILON_PLUS_DICHOTOM:
                return epsilon + sigma_n;
            default:
                return 1.0;
        }
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_LANGEVINSOLVER_H
