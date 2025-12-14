#ifndef DICHOTOMIC_NOISE_GENERATOR_H
#define DICHOTOMIC_NOISE_GENERATOR_H

#pragma once

#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <iostream>

namespace special {

/**
 * @class DichotomousNoiseGenerator
 * @brief Генератор дихотомного шума
 *
 * Реализует алгоритм из части 1 лабораторной работы
 * для создания коррелированного дихотомного процесса σ(t) ∈ {a, b}
 *
 * Свойства:
 * - ⟨σ(t)⟩ = 0 (нулевое среднее)
 * - ⟨σ(t)σ(t')⟩ = Q·exp(-Γ|t-t'|) (экспоненциальная корреляция)
 */
class DichotomousNoiseGenerator {
public:
    /**
     * @brief Конструктор
     * @param a, b значения шума (обычно +1 и -1)
     * @param gamma_a частота переходов a→b
     * @param gamma_b частота переходов b→a
     * @param seed для воспроизводимости
     */
    DichotomousNoiseGenerator(double a = 1.0,
                             double b = -1.0,
                             double gamma_a = 1.0,
                             double gamma_b = 1.0,
                             unsigned int seed = 42)
        : a_(a), b_(b),
          gamma_a_(gamma_a), gamma_b_(gamma_b),
          Gamma_(gamma_a + gamma_b),
          rng_(seed),
          uniform_dist_(0.0, 1.0)
    {
        if (Gamma_ <= 0.0) {
            throw std::invalid_argument("Gamma = gamma_a + gamma_b must be > 0");
        }
        // Проверка условия нулевого среднего
        if (std::abs(a + b) > 1e-6) {
            std::cerr << "Warning: For zero mean, require a + b ≈ 0\n"
                     << "Current: a=" << a << ", b=" << b
                     << ", a+b=" << (a+b) << "\n";
        }
        tau_c_ = 1.0 / Gamma_;
    }

    /**
     * @brief Генерировать реализацию дихотомного шума
     * @param dt временной шаг (должен быть << tau_c_)
     * @param n_steps количество шагов
     * @return вектор значений σ(t_i)
     */
    std::vector<double> generate(double dt, int n_steps)
    {
        if (dt <= 0.0) {
            throw std::invalid_argument("dt must be > 0");
        }
        if (n_steps < 1) {
            throw std::invalid_argument("n_steps must be >= 1");
        }

        // Проверка условия τ_C >> Δt
        if (dt > 0.1 * tau_c_) {
            std::cerr << "Warning: dt=" << dt << " may be too large compared to tau_c="
                     << tau_c_ << "\n"
                     << "Recommend: dt < 0.1 * tau_c = " << (0.1 * tau_c_) << "\n";
        }

        std::vector<double> sigma(n_steps);

        // Начальное значение (выбрать случайно или a)
        // По условию стационарности: P(a) = γ_b/Γ, P(b) = γ_a/Γ
        double prob_a = gamma_b_ / Gamma_;
        double sigma_current = (uniform_dist_(rng_) < prob_a) ? a_ : b_;
        sigma[0] = sigma_current;

        // Вычислить матрицу переходов
        const double exp_minus_Gamma_dt = std::exp(-Gamma_ * dt);

        // P(a_{n+1} | a_n, Δt) = γ_b/Γ + (γ_a/Γ)·exp(-Γ·Δt)
        const double P_aa = gamma_b_ / Gamma_ + (gamma_a_ / Gamma_) * exp_minus_Gamma_dt;

        // P(b_{n+1} | a_n, Δt) = (1 - exp(-Γ·Δt))·γ_b/Γ
        const double P_ab = (1.0 - exp_minus_Gamma_dt) * gamma_b_ / Gamma_;

        // P(a_{n+1} | b_n, Δt) = (1 - exp(-Γ·Δt))·γ_a/Γ
        const double P_ba = (1.0 - exp_minus_Gamma_dt) * gamma_a_ / Gamma_;

        // P(b_{n+1} | b_n, Δt) = γ_a/Γ + (γ_b/Γ)·exp(-Γ·Δt)
        const double P_bb = gamma_a_ / Gamma_ + (gamma_b_ / Gamma_) * exp_minus_Gamma_dt;

        // Генерировать траекторию
        for (int i = 1; i < n_steps; ++i) {
            double r = uniform_dist_(rng_);

            if (std::abs(sigma_current - a_) < 1e-10) {
                // В состоянии a: остаться в a с вероятностью P_aa, перейти в b с P_ab
                sigma_current = (r < P_aa) ? a_ : b_;
            } else {
                // В состоянии b: остаться в b с вероятностью P_bb, перейти в a с P_ba
                sigma_current = (r < P_ba) ? a_ : b_;
            }

            sigma[i] = sigma_current;
        }

        return sigma;
    }

    /**
     * @brief Вычислить корреляционную функцию численно
     * @param sigma вектор значений σ(t)
     * @return нормированная корреляция C(τ) = ⟨σ(t)σ(t+τ)⟩ / ⟨σ²⟩
     */
    std::vector<double> correlation(const std::vector<double>& sigma)
    {
        if (sigma.empty()) {
            throw std::invalid_argument("sigma vector is empty");
        }

        int n = static_cast<int>(sigma.size());
        std::vector<double> corr(n, 0.0);

        // Вычислить среднее σ (должно быть ≈ 0)
        double mean_sigma = std::accumulate(sigma.begin(), sigma.end(), 0.0) / n;

        // Вычислить дисперсию ⟨σ²⟩ - ⟨σ⟩²
        double var_sigma = 0.0;
        for (double s : sigma) {
            var_sigma += (s - mean_sigma) * (s - mean_sigma);
        }
        var_sigma /= n;

        if (var_sigma < 1e-15) {
            throw std::runtime_error("sigma has zero variance");
        }

        // Корреляция: C(τ) = ⟨σ(t)σ(t+τ)⟩ / var
        for (int tau = 0; tau < n; ++tau) {
            double sum = 0.0;
            int count = 0;
            for (int t = 0; t < n - tau; ++t) {
                sum += (sigma[t] - mean_sigma) * (sigma[t + tau] - mean_sigma);
                count++;
            }
            if (count > 0) {
                corr[tau] = sum / count / var_sigma;
            }
        }

        return corr;
    }

    /**
     * @brief Вычислить первый момент ⟨σ⟩
     */
    double mean(const std::vector<double>& sigma) const {
        if (sigma.empty()) return 0.0;
        return std::accumulate(sigma.begin(), sigma.end(), 0.0) / sigma.size();
    }

    /**
     * @brief Вычислить второй момент ⟨σ²⟩
     */
    double second_moment(const std::vector<double>& sigma) const {
        if (sigma.empty()) return 0.0;
        double sum = 0.0;
        for (double s : sigma) {
            sum += s * s;
        }
        return sum / sigma.size();
    }

    // ===== ПАРАМЕТРЫ ДОСТУПА =====

    double tau_c() const { return tau_c_; }
    double Gamma() const { return Gamma_; }
    double a() const { return a_; }
    double b() const { return b_; }
    double gamma_a() const { return gamma_a_; }
    double gamma_b() const { return gamma_b_; }

    /**
     * @brief Интенсивность шума Q
     */
    double Q() const {
        return a_ * b_ * Gamma_;
    }

    /**
     * @brief Проверить условия согласованности
     */
    void validate() const {
        std::cout << "DichotomousNoiseGenerator validation:\n";
        std::cout << "  a = " << a_ << ", b = " << b_ << "\n";
        std::cout << "  γ_a = " << gamma_a_ << ", γ_b = " << gamma_b_ << "\n";
        std::cout << "  Γ = γ_a + γ_b = " << Gamma_ << "\n";
        std::cout << "  τ_C = 1/Γ = " << tau_c_ << "\n";
        std::cout << "  Q = ab·Γ = " << Q() << "\n";
        std::cout << "  a + b = " << (a_ + b_) << " (should be ≈ 0)\n";
        std::cout << "  Asymmetry θ = |a - b| = " << std::abs(a_ - b_) << "\n";
    }

private:
    double a_, b_;
    double gamma_a_, gamma_b_;
    double Gamma_;  // = gamma_a + gamma_b = 1/tau_c
    double tau_c_;  // = 1 / Gamma
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_dist_;
};

} // namespace special

#endif // DICHOTOMIC_NOISE_GENERATOR_H