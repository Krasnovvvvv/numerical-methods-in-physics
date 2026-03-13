#ifndef NUMERICAL_METHODS_IN_PHYSICS_DIMLESSLANGEVINSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_DIMLESSLANGEVINSOLVER_H

#pragma once

#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"
#include "Labs/Special/Additional/Base/RatchetParams.h"

struct EnsembleResult {
    std::vector<double> sim_time;
    std::size_t n_particles = 0;
    bool has_trajectory = false;

    std::vector<double> mean_x;
    double mean_velocity = 0.0;
    double mean_x_final = 0.0;
};

class CyclicBarrier {
public:
    explicit CyclicBarrier(std::size_t count)
        : threshold_(count), count_(count), generation_(0) {}

    void wait() {
        std::unique_lock<std::mutex> lock(mutex_);
        const std::size_t gen = generation_;

        if (--count_ == 0) {
            generation_++;
            count_ = threshold_;
            cv_.notify_all();
            return;
        }

        cv_.wait(lock, [this, gen]() { return generation_ != gen; });
    }

private:
    std::mutex mutex_;
    std::condition_variable cv_;
    const std::size_t threshold_;
    std::size_t count_;
    std::size_t generation_;
};

struct alignas(64) PaddedDouble {
    double value = 0.0;
};

class DimLessLangevinSolver {
public:
    explicit DimLessLangevinSolver(const RatchetParams& params,
                                   unsigned int seed = std::random_device{}())
        : params_(params),
          rng_(seed),
          dt_(params.dt),
          sqrt_2dt_(std::sqrt(2.0 * params.dt))
    {
        if (params_.dt <= 0.0) {
            throw std::invalid_argument("dt must be > 0");
        }

        n_threads_ = std::thread::hardware_concurrency();
        if (n_threads_ == 0) n_threads_ = 4;
    }

    EnsembleResult solve(std::vector<DichotomicNoise>& noises,
                         std::size_t n_particles,
                         std::size_t total_time,
                         std::size_t burn_in = 0,
                         double x0 = 0.0,
                         bool store_trajectory = true,
                         std::size_t trajectory_stride = 1)
    {
        if (n_particles == 0) {
            throw std::invalid_argument("n_particles must be > 0");
        }
        if (noises.size() < n_particles) {
            throw std::invalid_argument("noises.size() must be >= n_particles");
        }

        EnsembleResult result;
        result.n_particles = n_particles;
        result.has_trajectory = store_trajectory;

        const std::size_t N = static_cast<std::size_t>(total_time / params_.dt);
        if (N == 0) {
            return result;
        }

        if (burn_in == 0) burn_in = N / 10;
        if (trajectory_stride == 0) trajectory_stride = 1;

        std::vector<double> xu(n_particles, x0);
        std::vector<double> xw(n_particles, wrap_unit_(x0));

        if (store_trajectory) {
            const std::size_t reserve_size = N / trajectory_stride + 2;
            result.sim_time.reserve(reserve_size);
            result.mean_x.reserve(reserve_size);
        }

        const std::size_t n_active = std::min(n_threads_, n_particles);
        const std::size_t chunk = (n_particles + n_active - 1) / n_active;

        std::vector<PaddedDouble> partial_sums(n_active);
        std::vector<std::mt19937> thread_rngs(n_active);

        for (std::size_t t = 0; t < n_active; ++t) {
            thread_rngs[t].seed(rng_() + static_cast<unsigned int>(t * 7919u + 17u));
        }

        CyclicBarrier barrier(n_active + 1);
        std::vector<std::thread> workers;
        workers.reserve(n_active);

        for (std::size_t t = 0; t < n_active; ++t) {
            const std::size_t p_begin = t * chunk;
            const std::size_t p_end = std::min(n_particles, p_begin + chunk);

            workers.emplace_back([this,
                                  &xu,
                                  &xw,
                                  &noises,
                                  &partial_sums,
                                  &thread_rngs,
                                  &barrier,
                                  p_begin,
                                  p_end,
                                  t,
                                  N]() {
                std::normal_distribution<double> gauss(0.0, 1.0);
                auto& local_rng = thread_rngs[t];

                for (std::size_t n = 0; n < N; ++n) {
                    double local_sum = 0.0;

                    for (std::size_t p = p_begin; p < p_end; ++p) {
                        const double sigma = noises[p].step();
                        const double f_t = compute_f_(sigma);

                        const double w_n = gauss(local_rng);
                        const double noise = sqrt_2dt_ * w_n;

                        const double dUdx_n = f_t * params_.dVdx(xw[p]);
                        const double xu_tilde = xu[p] - dUdx_n * dt_ + noise;
                        const double dUdx_tilde = f_t * params_.dVdx(wrap_unit_(xu_tilde));

                        const double delta =
                            -0.5 * (dUdx_n + dUdx_tilde) * dt_ + noise;

                        xu[p] += delta;
                        xw[p] = wrap_unit_(xu[p]);

                        local_sum += xu[p];
                    }

                    partial_sums[t].value = local_sum;

                    barrier.wait();
                    barrier.wait();
                }
            });
        }

        long double sum_t = 0.0L;
        long double sum_x = 0.0L;
        long double sum_tt = 0.0L;
        long double sum_tx = 0.0L;
        std::size_t fit_count = 0;

        for (std::size_t n = 0; n < N; ++n) {
            barrier.wait();

            double total_sum = 0.0;
            for (std::size_t t = 0; t < n_active; ++t) {
                total_sum += partial_sums[t].value;
            }

            const double mean_x = total_sum / static_cast<double>(n_particles);
            const double time = static_cast<double>(n) * dt_;

            result.mean_x_final = mean_x;

            if (n >= burn_in) {
                const long double tl = static_cast<long double>(time);
                const long double xl = static_cast<long double>(mean_x);
                sum_t += tl;
                sum_x += xl;
                sum_tt += tl * tl;
                sum_tx += tl * xl;
                ++fit_count;
            }

            if (store_trajectory && ((n % trajectory_stride) == 0 || n + 1 == N)) {
                result.sim_time.push_back(time);
                result.mean_x.push_back(mean_x);
            }

            barrier.wait();
        }

        for (auto& th : workers) {
            if (th.joinable()) th.join();
        }

        if (fit_count >= 2) {
            const long double denom =
                static_cast<long double>(fit_count) * sum_tt - sum_t * sum_t;

            if (std::abs(static_cast<double>(denom)) > 0.0) {
                result.mean_velocity = static_cast<double>(
                    (static_cast<long double>(fit_count) * sum_tx - sum_t * sum_x) / denom
                );
            }
        }

        return result;
    }

private:
    static double compute_f_(const double sigma, const double alpha) noexcept {
        return (sigma >= 0.0) ? 1.0 : alpha;
    }

    double compute_f_(const double sigma) const noexcept {
        return compute_f_(sigma, params_.alpha);
    }

    static double wrap_unit_(double x) noexcept {
        x -= std::floor(x);
        return x;
    }

    RatchetParams params_;
    std::mt19937 rng_;
    std::size_t n_threads_ = 4;
    double dt_ = 0.0;
    double sqrt_2dt_ = 0.0;
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_DIMLESSLANGEVINSOLVER_H
