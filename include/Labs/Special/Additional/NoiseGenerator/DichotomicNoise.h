#ifndef NUMERICAL_METHODS_IN_PHYSICS_DICHOTOMICNOISE_H
#define NUMERICAL_METHODS_IN_PHYSICS_DICHOTOMICNOISE_H

#pragma once

#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

struct DichotomicProfile {
    std::vector<double> t;
    std::vector<double> s;
    double m1;
    double m2;
    double m1_theory;
    double tau_c;
};

class DichotomicNoise {
public:
    enum class Mode {
        Symmetric,
        General,
        ZeroMean
    };

    /**
     * Конструктор для СИММЕТРИЧНОГО нулесреднего шума
     */
    DichotomicNoise(double a,
                    double tau_c,
                    double dt,
                    unsigned int seed = std::random_device{}())
        : a_(a),
          b_(-a),
          tau_c_(tau_c),
          dt_(dt),
          mode_(Mode::Symmetric),
          rng_(seed),
          uni_(0.0, 1.0)
    {
        if (tau_c_ <= 0.0) throw std::invalid_argument("tau_c must be > 0");
        if (dt_ <= 0.0) throw std::invalid_argument("dt must be > 0");
        if (a_ == 0.0) throw std::invalid_argument("a must be != 0");
        const double Gamma = 1.0 / tau_c_;
        gamma_a_ = 0.5 * Gamma;
        gamma_b_ = 0.5 * Gamma;
        init_transition_probs();
        init_state();
    }

    /**
     * Конструктор для АСИММЕТРИЧНОГО шума с ГАРАНТИРОВАННЫМ нулевым средним
     */
    static DichotomicNoise ZeroMeanAsymmetric(
        double a,
        double b,
        double gamma_b,
        double dt,
        unsigned int seed = std::random_device{}())
    {
        if (std::abs(b) < 1e-14) {
            throw std::invalid_argument("b must be != 0 for zero-mean asymmetric noise");
        }

        double gamma_a = -(a / b) * gamma_b;
        if (gamma_a <= 0.0) {
            throw std::invalid_argument(
                "With given a and b, computed gamma_a < 0. "
                "Try: a and b with opposite signs.");
        }

        DichotomicNoise gen(a, b, gamma_a, gamma_b, dt, seed);
        gen.mode_ = Mode::ZeroMean;
        return gen;
    }

    /**
     * Конструктор для ОБЩЕГО асимметричного шума
     */
    DichotomicNoise(double a,
                    double b,
                    double gamma_a,
                    double gamma_b,
                    double dt,
                    unsigned int seed = std::random_device{}())
        : a_(a),
          b_(b),
          gamma_a_(gamma_a),
          gamma_b_(gamma_b),
          tau_c_(1.0 / (gamma_a + gamma_b)),
          dt_(dt),
          mode_(Mode::General),
          rng_(seed),
          uni_(0.0, 1.0)
    {
        if (gamma_a_ <= 0.0 || gamma_b_ <= 0.0)
            throw std::invalid_argument("gamma_a and gamma_b must be > 0");
        if (dt_ <= 0.0)
            throw std::invalid_argument("dt must be > 0");
        init_transition_probs();
        init_state();
    }

    /**
     * Один шаг интегрирования: обновляет внутреннее состояние и возвращает текущее значение σ
     */
    double step() {
        double r = uni_(rng_);
        if (sigma_current_ == a_) {
            if (r > P_aa_) {
                sigma_current_ = b_;
            }
        } else {
            if (r > P_bb_) {
                sigma_current_ = a_;
            }
        }
        return sigma_current_;
    }

    /**
     * Сгенерировать полный профиль (для анализа статистики)
     */
    DichotomicProfile generate(std::size_t N,
                               std::size_t burn_in = 0,
                               bool start_from_stationary = true)
    {
        DichotomicProfile res;
        res.t.reserve(N);
        res.s.reserve(N);
        res.tau_c = tau_c_;

        // Выбор начального состояния
        double sigma;
        if (start_from_stationary) {
            double Gamma = gamma_a_ + gamma_b_;
            double P_a_stat = gamma_b_ / Gamma;
            double r0 = uni_(rng_);
            sigma = (r0 < P_a_stat) ? a_ : b_;
        } else {
            sigma = a_;
        }

        // Генерирование траектории
        for (std::size_t n = 0; n < N; ++n) {
            res.t.push_back(static_cast<double>(n) * dt_);
            res.s.push_back(sigma);

            // Марковский переход
            double r = uni_(rng_);
            if (sigma == a_) {
                if (r > P_aa_) {
                    sigma = b_;
                }
            } else {
                if (r > P_bb_) {
                    sigma = a_;
                }
            }
        }

        // Статистика после burn_in
        if (burn_in > N) burn_in = N;
        double sum1 = 0.0, sum2 = 0.0;
        std::size_t cnt = 0;
        for (std::size_t n = burn_in; n < N; ++n) {
            double v = res.s[n];
            sum1 += v;
            sum2 += v * v;
            ++cnt;
        }

        if (cnt > 0) {
            const double inv = 1.0 / static_cast<double>(cnt);
            res.m1 = sum1 * inv;
            res.m2 = sum2 * inv;
        } else {
            res.m1 = 0.0;
            res.m2 = 0.0;
        }

        double Gamma = gamma_a_ + gamma_b_;
        res.m1_theory = (a_ * gamma_b_ + b_ * gamma_a_) / Gamma;

        return res;
    }

    std::vector<double> autocorr_normalized(const std::vector<double>& s,
                                            std::size_t max_lag) const
    {
        // предполагаем ~ 0
        std::size_t N = s.size();
        if (max_lag >= N) max_lag = N - 1;

        double m2 = 0.0;
        for (double v : s) m2 += v * v;
        m2 /= static_cast<double>(N);
        if (m2 == 0.0) return std::vector<double>(max_lag + 1, 0.0);

        std::vector<double> C(max_lag + 1);
        C[0] = 1.0;

        for (std::size_t k = 1; k <= max_lag; ++k) {
            double sum = 0.0;
            std::size_t cnt = N - k;
            for (std::size_t n = 0; n < cnt; ++n) {
                sum += s[n] * s[n + k];
            }

            double ck = sum / static_cast<double>(cnt);
            C[k] = ck / m2;
        }

        return C;
    }

    double a() const { return a_; }
    double b() const { return b_; }
    double gamma_a() const { return gamma_a_; }
    double gamma_b() const { return gamma_b_; }
    double tau_c() const { return tau_c_; }

    double mean_theoretical() const {
        double Gamma = gamma_a_ + gamma_b_;
        return (a_ * gamma_b_ + b_ * gamma_a_) / Gamma;
    }

    double variance_theoretical() const {
        double Gamma = gamma_a_ + gamma_b_;
        double m1_th = (a_ * gamma_b_ + b_ * gamma_a_) / Gamma;
        double m2_th = (a_ * a_ * gamma_b_ + b_ * b_ * gamma_a_) / Gamma;
        return m2_th - m1_th * m1_th;
    }

    Mode mode() const { return mode_; }

private:
    double a_, b_;
    double gamma_a_, gamma_b_;
    double tau_c_, dt_;
    Mode mode_;
    double P_aa_; // P(a->a) = exp(-gamma_a * dt)
    double P_bb_; // P(b->b) = exp(-gamma_b * dt)
    double sigma_current_; // текущее состояние для on-the-fly генерации
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uni_;

    void init_transition_probs() {
        P_aa_ = std::exp(-gamma_a_ * dt_);
        P_bb_ = std::exp(-gamma_b_ * dt_);
    }

    void init_state() {
        double Gamma = gamma_a_ + gamma_b_;
        double P_a_stat = gamma_b_ / Gamma;
        double r0 = uni_(rng_);
        sigma_current_ = (r0 < P_a_stat) ? a_ : b_;
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_DICHOTOMICNOISE_H