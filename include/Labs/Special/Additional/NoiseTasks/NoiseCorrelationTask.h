#ifndef NUMERICAL_METHODS_IN_PHYSICS_NOISECORRELATIONTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_NOISECORRELATIONTASK_H

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "Labs/Special/Additional/NoiseTasks/INoiseTask.h"
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"
#include "Helpers/Plotter.h"

struct NoiseCorrelationResult {
    std::vector<double> t;
    std::vector<double> acf_exp;
    std::vector<double> acf_theory;

    double tau_c_theory = 0.0;
    double max_abs_error = 0.0;
    double rmse = 0.0;

    std::size_t n_trajectories = 0;
    std::size_t burn_in = 0;
    std::size_t measure_steps = 0;
    std::size_t max_lag = 0;
};

class NoiseCorrelationTask : public INoiseTask {
    double a_;
    double dt_;
    std::vector<double> tau_cs_;

    std::size_t burn_in_;
    std::size_t measure_steps_;
    std::size_t n_trajectories_;
    std::size_t max_lag_;

    unsigned int base_seed_;

    Plotter& plotter_;

    std::vector<DichotomicProfile> profiles_;
    std::vector<NoiseCorrelationResult> results_;

public:
    NoiseCorrelationTask(
        double a,
        double dt,
        const std::vector<double>& tau_cs,
        std::size_t burn_in,
        std::size_t measure_steps,
        std::size_t n_trajectories,
        std::size_t max_lag,
        Plotter& plotter,
        unsigned int base_seed = 700u)
        : a_(a),
          dt_(dt),
          tau_cs_(tau_cs),
          burn_in_(burn_in),
          measure_steps_(measure_steps),
          n_trajectories_(n_trajectories),
          max_lag_(max_lag),
          base_seed_(base_seed),
          plotter_(plotter)
    {}

    std::string name() const override {
        return "Calculating autocorrelation function for several correlation times";
    }

    void run() override {
        results_.clear();
        profiles_.clear();

        if (tau_cs_.empty() || measure_steps_ == 0 || n_trajectories_ == 0) {
            return;
        }

        results_.reserve(tau_cs_.size());
        profiles_.reserve(tau_cs_.size());

        for (std::size_t i = 0; i < tau_cs_.size(); ++i) {
            const double tau_c = tau_cs_[i];
            const unsigned int seed_shift =
                base_seed_ + static_cast<unsigned int>(10000u * i);

            results_.push_back(calculate_for_tau(tau_c, seed_shift));
        }

        plot_results();
        print();
    }

    const std::vector<NoiseCorrelationResult>& results() const {
        return results_;
    }

    const std::vector<DichotomicProfile>& profiles() const {
        return profiles_;
    }

    void print(std::size_t n_show = 10) const {
        std::cout << "=== Noise correlation (ensemble, several tau_c) ===\n";
        std::cout << "a                : " << a_ << '\n';
        std::cout << "dt               : " << dt_ << '\n';
        std::cout << "burn_in          : " << burn_in_ << '\n';
        std::cout << "measure_steps    : " << measure_steps_ << '\n';
        std::cout << "n_trajectories   : " << n_trajectories_ << '\n';
        std::cout << "max_lag          : " << max_lag_ << '\n';
        std::cout << "n_tau_values     : " << tau_cs_.size() << "\n\n";

        for (std::size_t i = 0; i < results_.size(); ++i) {
            const auto& r = results_[i];

            std::cout << "--- tau_c = " << r.tau_c_theory << " ---\n";
            std::cout << "max_abs_error    : " << r.max_abs_error << '\n';
            std::cout << "rmse             : " << r.rmse << '\n';

            const std::size_t m = std::min(n_show, r.acf_exp.size());
            std::cout << "lag\t"
                      << "t\t"
                      << "acf_exp\t"
                      << "acf_theory\n";

            for (std::size_t k = 0; k < m; ++k) {
                std::cout << k << '\t'
                          << r.t[k] << '\t'
                          << r.acf_exp[k] << '\t'
                          << r.acf_theory[k] << '\n';
            }

            std::cout << '\n';
        }
    }

private:
    NoiseCorrelationResult calculate_for_tau(double tau_c, unsigned int seed_base) {
        NoiseCorrelationResult result;
        DichotomicProfile profile;

        result.n_trajectories = n_trajectories_;
        result.burn_in = burn_in_;
        result.measure_steps = measure_steps_;

        if (measure_steps_ == 0 || n_trajectories_ == 0) {
            profiles_.push_back(profile);
            return result;
        }

        const std::size_t lag_eff = std::min(max_lag_, measure_steps_ - 1);
        result.max_lag = lag_eff;
        result.t.assign(lag_eff + 1, 0.0);
        result.acf_exp.assign(lag_eff + 1, 0.0);
        result.acf_theory.assign(lag_eff + 1, 0.0);

        for (std::size_t k = 0; k <= lag_eff; ++k) {
            result.t[k] = static_cast<double>(k) * dt_;
        }

        for (std::size_t tr = 0; tr < n_trajectories_; ++tr) {
            DichotomicNoise noise(
                a_,
                tau_c,
                dt_,
                seed_base + static_cast<unsigned int>(tr)
            );

            if (tr == 0) {
                result.tau_c_theory = noise.tau_c();

                profile.tau_c = noise.tau_c();
                profile.t.reserve(measure_steps_);
                profile.s.reserve(measure_steps_);
            }

            for (std::size_t i = 0; i < burn_in_; ++i) {
                noise.step();
            }

            std::vector<double> s;
            s.reserve(measure_steps_);

            for (std::size_t i = 0; i < measure_steps_; ++i) {
                const double sigma = noise.step();
                s.push_back(sigma);

                if (tr == 0) {
                    profile.t.push_back(static_cast<double>(i) * dt_);
                    profile.s.push_back(sigma);
                }
            }

            if (tr == 0 && !s.empty()) {
                long double sum1 = 0.0L;
                long double sum2 = 0.0L;

                for (double v : s) {
                    sum1 += static_cast<long double>(v);
                    sum2 += static_cast<long double>(v) * static_cast<long double>(v);
                }

                const long double inv = 1.0L / static_cast<long double>(s.size());
                profile.m1 = static_cast<double>(sum1 * inv);
                profile.m2 = static_cast<double>(sum2 * inv);
                profile.m1_theory = noise.mean_theoretical();
            }

            const std::vector<double> acf = noise.autocorr_normalized(s, lag_eff);

            for (std::size_t k = 0; k < acf.size(); ++k) {
                result.acf_exp[k] += acf[k];
            }
        }

        const double inv_tr = 1.0 / static_cast<double>(n_trajectories_);
        for (double& v : result.acf_exp) {
            v *= inv_tr;
        }

        long double sq_sum = 0.0L;
        for (std::size_t k = 0; k <= lag_eff; ++k) {
            result.acf_theory[k] = std::exp(-result.t[k] / result.tau_c_theory);

            const double err = std::abs(result.acf_exp[k] - result.acf_theory[k]);
            result.max_abs_error = std::max(result.max_abs_error, err);
            sq_sum += static_cast<long double>(err) * static_cast<long double>(err);
        }

        result.rmse = std::sqrt(
            static_cast<double>(sq_sum / static_cast<long double>(lag_eff + 1))
        );

        profiles_.push_back(profile);
        return result;
    }

    void plot_results() {
        if (results_.empty()) {
            return;
        }

        std::vector<std::string> labels;
        std::vector<std::vector<double>> xs;
        std::vector<std::vector<double>> ys;

        labels.reserve(results_.size() * 2);
        xs.reserve(results_.size() * 2);
        ys.reserve(results_.size() * 2);

        for (const auto& r : results_) {
            if (r.t.empty() || r.acf_exp.empty() || r.acf_theory.empty()) {
                continue;
            }

            labels.emplace_back(
                "Expected ACF, tau_c = " + std::to_string(r.tau_c_theory)
            );
            xs.emplace_back(r.t);
            ys.emplace_back(r.acf_exp);

            labels.emplace_back(
                "Theory exp(-t/tau_c), tau_c = " + std::to_string(r.tau_c_theory)
            );
            xs.emplace_back(r.t);
            ys.emplace_back(r.acf_theory);
        }

        if (!xs.empty()) {
            plotter_.plot(
                xs,
                ys,
                labels,
                "t",
                "C(t)"
            );
        }
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_NOISECORRELATIONTASK_H