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
    double m1 = 0.0;
    double m2 = 0.0;
    double m1_theory = 0.0;
    double tau_c = 0.0;
};

class DichotomicNoise {
public:
    enum class Mode {
        Symmetric,
        General,
        ZeroMean
    };

    // Симметричный нулесредний шум: a, -a, tau_c = 1 / Gamma
    DichotomicNoise(double a,
                    double tau_c,
                    double dt,
                    unsigned int seed = std::random_device{}())
        : a_(a),
          b_(-a),
          gamma_a_(0.0),
          gamma_b_(0.0),
          tau_c_(tau_c),
          dt_(dt),
          mode_(Mode::Symmetric),
          rng_(seed),
          uni_(0.0, 1.0)
    {
        if (a_ == 0.0) throw std::invalid_argument("a must be != 0");
        if (tau_c_ <= 0.0) throw std::invalid_argument("tau_c must be > 0");
        if (dt_ <= 0.0) throw std::invalid_argument("dt must be > 0");

        const double Gamma = 1.0 / tau_c_;
        gamma_a_ = 0.5 * Gamma;
        gamma_b_ = 0.5 * Gamma;

        init_transition_probs();
        init_state();
    }

    // Общий асимметричный шум
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
          tau_c_(0.0),
          dt_(dt),
          mode_(Mode::General),
          rng_(seed),
          uni_(0.0, 1.0)
    {
        if (gamma_a_ <= 0.0 || gamma_b_ <= 0.0) {
            throw std::invalid_argument("gamma_a and gamma_b must be > 0");
        }
        if (dt_ <= 0.0) {
            throw std::invalid_argument("dt must be > 0");
        }

        tau_c_ = 1.0 / (gamma_a_ + gamma_b_);

        init_transition_probs();
        init_state();
    }

    // Асимметричный шум с нулевым средним: a * gamma_b + b * gamma_a = 0
    static DichotomicNoise ZeroMeanAsymmetric(double a,
                                              double b,
                                              double gamma_b,
                                              double dt,
                                              unsigned int seed = std::random_device{}())
    {
        if (std::abs(b) < 1e-14) {
            throw std::invalid_argument("b must be != 0 for zero-mean asymmetric noise");
        }
        if (gamma_b <= 0.0) {
            throw std::invalid_argument("gamma_b must be > 0");
        }

        const double gamma_a = -(a / b) * gamma_b;
        if (gamma_a <= 0.0) {
            throw std::invalid_argument(
                "For zero-mean asymmetric noise, a and b must have opposite signs");
        }

        DichotomicNoise gen(a, b, gamma_a, gamma_b, dt, seed);
        gen.mode_ = Mode::ZeroMean;
        return gen;
    }

    double step() {
        const double r = uni_(rng_);

        if (state_is_a_) {
            if (r >= P_aa_) state_is_a_ = false;
        } else {
            if (r >= P_bb_) state_is_a_ = true;
        }

        return state_is_a_ ? a_ : b_;
    }

    double current() const noexcept {
        return state_is_a_ ? a_ : b_;
    }

    DichotomicProfile generate(std::size_t N,
                               std::size_t burn_in = 0,
                               bool start_from_stationary = true)
    {
        DichotomicProfile res;
        res.t.reserve(N);
        res.s.reserve(N);
        res.tau_c = tau_c_;

        bool local_state_is_a = start_from_stationary ? sample_stationary_state_bool() : true;

        for (std::size_t n = 0; n < N; ++n) {
            const double r = uni_(rng_);
            if (local_state_is_a) {
                if (r >= P_aa_) local_state_is_a = false;
            } else {
                if (r >= P_bb_) local_state_is_a = true;
            }

            res.t.push_back(static_cast<double>(n) * dt_);
            res.s.push_back(local_state_is_a ? a_ : b_);
        }

        if (burn_in > N) burn_in = N;

        double sum1 = 0.0;
        double sum2 = 0.0;
        std::size_t cnt = 0;

        for (std::size_t n = burn_in; n < N; ++n) {
            const double v = res.s[n];
            sum1 += v;
            sum2 += v * v;
            ++cnt;
        }

        if (cnt > 0) {
            const double inv = 1.0 / static_cast<double>(cnt);
            res.m1 = sum1 * inv;
            res.m2 = sum2 * inv;
        }

        res.m1_theory = mean_theoretical();
        return res;
    }

    std::vector<double> autocorr_normalized(const std::vector<double>& s,
                                            std::size_t max_lag) const
    {
        const std::size_t N = s.size();
        if (N == 0) return {};
        if (max_lag >= N) max_lag = N - 1;

        double mean = 0.0;
        for (double v : s) mean += v;
        mean /= static_cast<double>(N);

        double var = 0.0;
        for (double v : s) {
            const double dv = v - mean;
            var += dv * dv;
        }
        var /= static_cast<double>(N);

        if (var <= 0.0) {
            return std::vector<double>(max_lag + 1, 0.0);
        }

        std::vector<double> C(max_lag + 1, 0.0);
        C[0] = 1.0;

        for (std::size_t k = 1; k <= max_lag; ++k) {
            double sum = 0.0;
            const std::size_t cnt = N - k;
            for (std::size_t n = 0; n < cnt; ++n) {
                sum += (s[n] - mean) * (s[n + k] - mean);
            }
            C[k] = (sum / static_cast<double>(cnt)) / var;
        }

        return C;
    }

    double a() const noexcept { return a_; }
    double b() const noexcept { return b_; }
    double gamma_a() const noexcept { return gamma_a_; }
    double gamma_b() const noexcept { return gamma_b_; }
    double tau_c() const noexcept { return tau_c_; }
    Mode mode() const noexcept { return mode_; }

    double mean_theoretical() const {
        const double Gamma = gamma_a_ + gamma_b_;
        return (a_ * gamma_b_ + b_ * gamma_a_) / Gamma;
    }

    double variance_theoretical() const {
        const double Gamma = gamma_a_ + gamma_b_;
        const double m1 = mean_theoretical();
        const double m2 = (a_ * a_ * gamma_b_ + b_ * b_ * gamma_a_) / Gamma;
        return m2 - m1 * m1;
    }

private:
    double a_ = 0.0;
    double b_ = 0.0;
    double gamma_a_ = 0.0;
    double gamma_b_ = 0.0;
    double tau_c_ = 0.0;
    double dt_ = 0.0;

    Mode mode_ = Mode::General;

    double P_aa_ = 0.0;
    double P_bb_ = 0.0;
    double P_a_stat_ = 0.5;

    bool state_is_a_ = true;

    std::mt19937 rng_;
    std::uniform_real_distribution<double> uni_;

    void init_transition_probs() {
        const double Gamma = gamma_a_ + gamma_b_;
        const double e = std::exp(-Gamma * dt_);

        P_aa_ = gamma_b_ / Gamma + gamma_a_ / Gamma * e;
        P_bb_ = gamma_a_ / Gamma + gamma_b_ / Gamma * e;
        P_a_stat_ = gamma_b_ / Gamma;
    }

    bool sample_stationary_state_bool() {
        return uni_(rng_) < P_a_stat_;
    }

    void init_state() {
        state_is_a_ = sample_stationary_state_bool();
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_DICHOTOMICNOISE_H