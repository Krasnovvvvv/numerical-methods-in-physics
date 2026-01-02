#ifndef NUMERICAL_METHODS_IN_PHYSICS_LANGEVINSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_LANGEVINSOLVER_H

#pragma once

#include <cmath>
#include <functional>
#include <random>
#include <stdexcept>
#include <vector>
#include <thread>
#include <mutex>
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"

// ============================================================================
// Структура для одной траектории (старая версия, для совместимости)
// ============================================================================
struct LangevinTrajectory {
    std::vector<double> t; // время
    std::vector<double> x; // позиция
    double L; // период
    double mean_position;
    double mean_velocity;
    double diffusion_coeff;
};

// ============================================================================
// Структура для ансамбля частиц
// ============================================================================
struct LangevinEnsembleResult {
    std::vector<double> t; // временные точки (хранится, если нужна траектория)
    std::vector<double> mean_x; // средняя позиция по ансамблю
    double L; // период
    double mean_velocity; // средняя скорость по ансамблю
    double diffusion_coeff; // коэффициент диффузии
    std::size_t n_particles; // число частиц
    bool has_trajectory; // флаг: хранится ли полная траектория
};

// ============================================================================
// Решатель Ланжевена
// ============================================================================
class LangevinSolver {
public:
    using DVDXFunc = std::function<double(double)>;

    enum class ModulationType {
        CONSTANT, // f(t) = 1.0
        DICHOTOM_SYMMETRIC, // f(t) = (1/2)(1 + σ(t))
        EPSILON_PLUS_DICHOTOM // f(t) = ε + σ(t)
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
          gaussian_(0.0, 1.0),
          sqrt2_coeff_(std::sqrt(2.0)),
          L_half_(L / 2.0),
          n_threads_(std::thread::hardware_concurrency())
    {
        if (n_threads_ == 0) n_threads_ = 4;
        if (dt_ <= 0.0) throw std::invalid_argument("dt must be > 0");
        if (L_ <= 0.0) throw std::invalid_argument("L must be > 0");
        if (V0_ < 0.0) throw std::invalid_argument("V0 must be >= 0");
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
        res.t.reserve(N);
        res.x.reserve(N);
        res.L = L_;

        double x = std::fmod(x0, L_);
        if (x < 0.0) x += L_;

        for (std::size_t n = 0; n < N; ++n) {
            res.t.push_back(static_cast<double>(n) * dt_);
            res.x.push_back(x);

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
                if (dx > L_half_) dx -= L_;
                if (dx < -L_half_) dx += L_;
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
                if (dx > L_half_) dx -= L_;
                if (dx < -L_half_) dx += L_;
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
    // Ансамбль с общим шумом (многопоточная версия)
    // ========================================================================
    LangevinEnsembleResult solve_ensemble(DichotomicNoise& noise,
                                          std::size_t N,
                                          std::size_t n_particles,
                                          double x0 = 0.0,
                                          double epsilon = 0.0,
                                          std::size_t burn_in = 0,
                                          bool store_trajectory = true)
    {
        LangevinEnsembleResult res;
        res.L = L_;
        res.n_particles = n_particles;
        res.has_trajectory = store_trajectory;

        if (store_trajectory) {
            res.t.reserve(N);
            res.mean_x.reserve(N);
        }

        // Создаём копию шума для каждого потока
        std::vector<DichotomicNoise> noise_copies;
        for (std::size_t t = 0; t < n_threads_; ++t) {
            noise_copies.push_back(noise);
        }

        // Массив позиций частиц
        std::vector<double> x(n_particles);
        for (auto& xi : x) {
            xi = std::fmod(x0, L_);
            if (xi < 0.0) xi += L_;
        }

        // Генераторы случайных чисел для каждого потока
        std::vector<std::mt19937> thread_rngs(n_threads_);
        for (std::size_t t = 0; t < n_threads_; ++t) {
            thread_rngs[t].seed(rng_() + t);
        }

        std::vector<std::normal_distribution<double>> thread_gaussians(n_threads_, std::normal_distribution<double>(0.0, 1.0));

        // Интегрирование с многопоточностью
        double sum_msd = 0.0;
        std::size_t cnt_msd = 0;
        double prev_mean_x = 0.0;

        for (std::size_t n = 0; n < N; ++n) {
            double sigma_n = noise_copies[0].step();
            double f_t = compute_f_t(sigma_n, epsilon);

            // Разбиваем частицы на чанки для потоков
            std::size_t chunk = (n_particles + n_threads_ - 1) / n_threads_;
            std::vector<std::thread> threads;
            std::vector<double> partial_sums(n_threads_, 0.0);
            std::mutex sum_mutex;

            for (std::size_t t = 0; t < n_threads_; ++t) {
                std::size_t p_begin = t * chunk;
                std::size_t p_end = std::min(n_particles, p_begin + chunk);
                if (p_begin >= p_end) break;

                threads.emplace_back([&, t, p_begin, p_end]() {
                    double local_sum = 0.0;
                    auto& local_rng = thread_rngs[t];
                    auto& local_gaussian = thread_gaussians[t];

                    for (std::size_t p = p_begin; p < p_end; ++p) {
                        double zeta_n = local_gaussian(local_rng);
                        double dVdx_n = dVdx_(x[p]);
                        double dxdt = -V0_ * dVdx_n * f_t + sqrt2_coeff_ * zeta_n;

                        x[p] += dxdt * dt_;
                        x[p] = std::fmod(x[p], L_);
                        if (x[p] < 0.0) x[p] += L_;

                        local_sum += x[p];
                    }

                    {
                        std::lock_guard<std::mutex> lock(sum_mutex);
                        partial_sums[t] = local_sum;
                    }
                });
            }

            // Ждём завершения всех потоков
            for (auto& th : threads) {
                if (th.joinable()) th.join();
            }

            // Агрегируем суммы
            double sum_x = 0.0;
            for (double s : partial_sums) sum_x += s;

            double mean_x_n = sum_x / static_cast<double>(n_particles);

            if (store_trajectory) {
                res.t.push_back(static_cast<double>(n) * dt_);
                res.mean_x.push_back(mean_x_n);
            }

            // Накопление MSD после burn_in для диффузии
            if (n > burn_in && n > 0) {
                double dmean_x = mean_x_n - prev_mean_x;
                if (dmean_x > L_half_) dmean_x -= L_;
                if (dmean_x < -L_half_) dmean_x += L_;
                sum_msd += dmean_x * dmean_x;
                ++cnt_msd;
            }

            if (n >= burn_in) {
                prev_mean_x = mean_x_n;
            }
        }

        // === Средняя скорость по ансамблю (после burn_in) ===
        double t_max = (N - burn_in) * dt_;
        if (N > burn_in + 1 && t_max > 0.0 && store_trajectory) {
            double sum_all_disp = 0.0;
            double x0_eff_mean = res.mean_x[0];
            for (std::size_t i = 0; i < res.mean_x.size(); ++i) {
                double dx = res.mean_x[i] - x0_eff_mean;
                if (dx > L_half_) dx -= L_;
                if (dx < -L_half_) dx += L_;
                sum_all_disp += dx;
            }

            std::size_t N_eff = N - burn_in;
            res.mean_velocity = (1.0 / t_max) * (1.0 / static_cast<double>(N_eff)) * sum_all_disp;
        } else {
            res.mean_velocity = 0.0;
        }

        // === Коэффициент диффузии ===
        if (cnt_msd > 0) {
            double mean_msd = sum_msd / static_cast<double>(cnt_msd);
            res.diffusion_coeff = mean_msd / (2.0 * dt_);
        } else {
            res.diffusion_coeff = 0.0;
        }

        return res;
    }

    // ========================================================================
    // Ансамбль с независимыми шумами (многопоточная версия)
    // ========================================================================
    LangevinEnsembleResult solve_ensemble_independent(
        std::vector<DichotomicNoise>& noises,
        std::size_t N,
        std::size_t n_particles,
        double x0 = 0.0,
        double epsilon = 0.0,
        std::size_t burn_in = 0,
        bool store_trajectory = true)
    {
        LangevinEnsembleResult res;
        res.L = L_;
        res.n_particles = n_particles;
        res.has_trajectory = store_trajectory;

        if (store_trajectory) {
            res.t.reserve(N);
            res.mean_x.reserve(N);
        }

        // Массив позиций частиц
        std::vector<double> x(n_particles);
        for (auto& xi : x) {
            xi = std::fmod(x0, L_);
            if (xi < 0.0) xi += L_;
        }

        // Генераторы случайных чисел для каждого потока
        std::vector<std::mt19937> thread_rngs(n_threads_);
        for (std::size_t t = 0; t < n_threads_; ++t) {
            thread_rngs[t].seed(rng_() + t);
        }

        std::vector<std::normal_distribution<double>> thread_gaussians(n_threads_, std::normal_distribution<double>(0.0, 1.0));

        // Интегрирование с многопоточностью
        double sum_msd = 0.0;
        std::size_t cnt_msd = 0;
        double prev_mean_x = 0.0;

        for (std::size_t n = 0; n < N; ++n) {
            // Разбиваем частицы на чанки для потоков
            std::size_t chunk = (n_particles + n_threads_ - 1) / n_threads_;
            std::vector<std::thread> threads;
            std::vector<double> partial_sums(n_threads_, 0.0);
            std::mutex sum_mutex;

            for (std::size_t t = 0; t < n_threads_; ++t) {
                std::size_t p_begin = t * chunk;
                std::size_t p_end = std::min(n_particles, p_begin + chunk);
                if (p_begin >= p_end) break;

                threads.emplace_back([&, t, p_begin, p_end]() {
                    double local_sum = 0.0;
                    auto& local_rng = thread_rngs[t];
                    auto& local_gaussian = thread_gaussians[t];

                    for (std::size_t p = p_begin; p < p_end; ++p) {
                        double sigma_n = noises[p].step(); // один шаг шума для частицы p
                        double f_t = compute_f_t(sigma_n, epsilon);
                        double zeta_n = local_gaussian(local_rng);
                        double dVdx_n = dVdx_(x[p]);
                        double dxdt = -V0_ * dVdx_n * f_t + sqrt2_coeff_ * zeta_n;

                        x[p] += dxdt * dt_;
                        x[p] = std::fmod(x[p], L_);
                        if (x[p] < 0.0) x[p] += L_;

                        local_sum += x[p];
                    }

                    {
                        std::lock_guard<std::mutex> lock(sum_mutex);
                        partial_sums[t] = local_sum;
                    }
                });
            }

            // Ждём завершения всех потоков
            for (auto& th : threads) {
                if (th.joinable()) th.join();
            }

            // Агрегируем суммы
            double sum_x = 0.0;
            for (double s : partial_sums) sum_x += s;

            double mean_x_n = sum_x / static_cast<double>(n_particles);

            if (store_trajectory) {
                res.t.push_back(static_cast<double>(n) * dt_);
                res.mean_x.push_back(mean_x_n);
            }

            // Накопление MSD после burn_in для диффузии
            if (n > burn_in && n > 0) {
                double dmean_x = mean_x_n - prev_mean_x;
                if (dmean_x > L_half_) dmean_x -= L_;
                if (dmean_x < -L_half_) dmean_x += L_;
                sum_msd += dmean_x * dmean_x;
                ++cnt_msd;
            }

            if (n >= burn_in) {
                prev_mean_x = mean_x_n;
            }
        }

        // === Средняя скорость (после burn_in) ===
        double t_max = (N - burn_in) * dt_;
        if (N > burn_in + 1 && t_max > 0.0 && store_trajectory) {
            double sum_all_disp = 0.0;
            double x0_eff_mean = res.mean_x[0];
            for (std::size_t i = 0; i < res.mean_x.size(); ++i) {
                double dx = res.mean_x[i] - x0_eff_mean;
                if (dx > L_half_) dx -= L_;
                if (dx < -L_half_) dx += L_;
                sum_all_disp += dx;
            }

            std::size_t N_eff = N - burn_in;
            res.mean_velocity = (1.0 / t_max) * (1.0 / static_cast<double>(N_eff)) * sum_all_disp;
        } else {
            res.mean_velocity = 0.0;
        }

        // === Коэффициент диффузии ===
        if (cnt_msd > 0) {
            double mean_msd = sum_msd / static_cast<double>(cnt_msd);
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
    double L_half_;
    std::mt19937 rng_;
    std::normal_distribution<double> gaussian_;
    std::size_t n_threads_;

    // Вспомогательная функция для вычисления модуляции f(t)
    inline double compute_f_t(double sigma_n, double epsilon) const {
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