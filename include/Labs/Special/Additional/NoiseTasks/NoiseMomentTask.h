#ifndef NUMERICAL_METHODS_IN_PHYSICS_NOISEMOMENTTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_NOISEMOMENTTASK_H

#pragma once

#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>

#include "Labs/Special/Additional/NoiseTasks/INoiseTask.h"
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"

struct NoiseMomentsResult {
    double m1_exp = 0.0;
    double m2_exp = 0.0;

    double m1_theory = 0.0;
    double m2_theory = 0.0;

    double var_exp = 0.0;
    double var_theory = 0.0;

    double abs_err_m1 = 0.0;
    double abs_err_m2 = 0.0;
    double abs_err_var = 0.0;

    std::size_t n_trajectories = 0;
    std::size_t burn_in = 0;
    std::size_t measure_steps = 0;
};

class NoiseMomentTask : public INoiseTask {
    double a_;
    double dt_;
    double tau_c_;

    std::size_t burn_in_;
    std::size_t measure_steps_;
    std::size_t n_trajectories_;
    unsigned int base_seed_;

    DichotomicProfile profile_;
    NoiseMomentsResult result_;

    void print() const {
        std::cout << "=== Noise moments (ensemble) ===\n";
        std::cout << "a                : " << a_ << '\n';
        std::cout << "tau_c            : " << tau_c_ << '\n';
        std::cout << "dt               : " << dt_ << '\n';
        std::cout << "burn_in          : " << burn_in_ << '\n';
        std::cout << "measure_steps    : " << measure_steps_ << '\n';
        std::cout << "n_trajectories   : " << n_trajectories_ << '\n';

        std::cout << "m1 exp           : " << result_.m1_exp << '\n';
        std::cout << "m1 theory        : " << result_.m1_theory << '\n';
        std::cout << "|m1 exp-theory|  : " << result_.abs_err_m1 << '\n';

        std::cout << "m2 exp           : " << result_.m2_exp << '\n';
        std::cout << "m2 theory        : " << result_.m2_theory << '\n';
        std::cout << "|m2 exp-theory|  : " << result_.abs_err_m2 << '\n';

        std::cout << "var exp          : " << result_.var_exp << '\n';
        std::cout << "var theory       : " << result_.var_theory << '\n';
        std::cout << "|var exp-theory| : " << result_.abs_err_var << '\n';
    }

public:
    NoiseMomentTask(
        double a,
        double dt,
        double tau_c,
        std::size_t burn_in,
        std::size_t measure_steps,
        std::size_t n_trajectories,
        unsigned int base_seed = 600u)
        : a_(a),
          dt_(dt),
          tau_c_(tau_c),
          burn_in_(burn_in),
          measure_steps_(measure_steps),
          n_trajectories_(n_trajectories),
          base_seed_(base_seed)
    {}

    std::string name() const override {
        return "Calculating the first and second moments";
    }

    void run() override {
        result_ = {};
        result_.n_trajectories = n_trajectories_;
        result_.burn_in = burn_in_;
        result_.measure_steps = measure_steps_;

        long double sum1 = 0.0L;
        long double sum2 = 0.0L;
        std::size_t count = 0;

        for (std::size_t tr = 0; tr < n_trajectories_; ++tr) {
            DichotomicNoise noise(
                a_,
                tau_c_,
                dt_,
                base_seed_ + static_cast<unsigned int>(tr)
            );

            if (tr == 0) {
                profile_ = noise.generate(burn_in_ + measure_steps_);
                result_.m1_theory = noise.mean_theoretical();
                result_.var_theory = noise.variance_theoretical();
                result_.m2_theory =
                    result_.var_theory + result_.m1_theory * result_.m1_theory;
            }

            for (std::size_t i = 0; i < burn_in_; ++i) {
                noise.step();
            }

            for (std::size_t i = 0; i < measure_steps_; ++i) {
                const double sigma = noise.step();
                sum1 += static_cast<long double>(sigma);
                sum2 += static_cast<long double>(sigma) * static_cast<long double>(sigma);
                ++count;
            }
        }

        if (count > 0) {
            const long double inv = 1.0L / static_cast<long double>(count);

            result_.m1_exp = static_cast<double>(sum1 * inv);
            result_.m2_exp = static_cast<double>(sum2 * inv);
            result_.var_exp = result_.m2_exp - result_.m1_exp * result_.m1_exp;

            result_.abs_err_m1 = std::abs(result_.m1_exp - result_.m1_theory);
            result_.abs_err_m2 = std::abs(result_.m2_exp - result_.m2_theory);
            result_.abs_err_var = std::abs(result_.var_exp - result_.var_theory);
        }
        print();
    }

    const NoiseMomentsResult& result() const {
        return result_;
    }

    const DichotomicProfile& profile() const {
        return profile_;
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_NOISEMOMENTTASK_H