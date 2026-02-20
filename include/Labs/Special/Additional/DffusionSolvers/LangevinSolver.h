#ifndef NUMERICAL_METHODS_IN_PHYSICS_LANGEVINSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_LANGEVINSOLVER_H

#pragma once

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <random>
#include <cmath>
#include <atomic>
#include <stdexcept>
#include <algorithm>
#include <memory>

#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"

struct LangevinEnsembleResult {
    std::vector<double> t;

    // Теперь mean_x = UNWRAPPED <x(t)> (для сравнения с теорией v*t)
    std::vector<double> mean_x;

    // Дополнительно: WRAPPED <x(t)> в [0, L) (удобно для контроля)
    std::vector<double> mean_x_wrapped;

    double L = 0.0;
    double mean_velocity = 0.0;
    double diffusion_coeff = 0.0;

    std::size_t n_particles = 0;
    bool has_trajectory = true;
};

class WorkerThread {
public:
    WorkerThread() : should_stop_(false), is_working_(false), is_running_(true) {
        thread_ = std::thread(&WorkerThread::run, this);
    }

    WorkerThread(const WorkerThread&) = delete;
    WorkerThread& operator=(const WorkerThread&) = delete;
    WorkerThread(WorkerThread&&) = delete;
    WorkerThread& operator=(WorkerThread&&) = delete;

    ~WorkerThread() { stop(); }

    void submit_task(std::function<void()> task) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!is_running_) return;
        task_queue_.push(std::move(task));
        cv_.notify_one();
    }

    void stop() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!is_running_) return;
        should_stop_ = true;
        is_running_ = false;
        cv_.notify_one();
        if (thread_.joinable()) thread_.join();
    }

private:
    void run() {
        while (true) {
            std::function<void()> task;

            {
                std::unique_lock<std::mutex> lock(mutex_);
                cv_.wait(lock, [this]() { return !task_queue_.empty() || should_stop_; });
                if (should_stop_ && task_queue_.empty()) break;

                task = std::move(task_queue_.front());
                task_queue_.pop();
                is_working_ = true;
            }

            if (task) task();

            {
                std::lock_guard<std::mutex> lock(mutex_);
                is_working_ = false;
                cv_.notify_one();
            }
        }
    }

    std::thread thread_;
    std::queue<std::function<void()>> task_queue_;
    std::mutex mutex_;
    std::condition_variable cv_;
    bool should_stop_;
    bool is_working_;
    bool is_running_;
};

class LangevinSolver {
public:
    using DVDXFunc = std::function<double(double)>;

    enum class ModulationType {
        CONSTANT,
        TWO_LEVEL_F
    };

    LangevinSolver(
        DVDXFunc dVdx,
        double L,
        double dt,
        double D,
        double beta,
        ModulationType mod_type = ModulationType::CONSTANT,
        unsigned int seed = std::random_device{}())
        : dVdx_(std::move(dVdx)),
          L_(L),
          dt_(dt),
          mod_type_(mod_type),
          D_(D),
          beta_(beta),
          rng_(seed),
          sqrt2D_coeff_(std::sqrt(2.0 * D_)),
          sqrt_dt_(std::sqrt(dt_)),
          L_half_(L_ / 2.0),
          n_threads_(std::thread::hardware_concurrency()) {
        if (n_threads_ == 0) n_threads_ = 4;
        validate_params();
        initialize_thread_pool();
    }

    ~LangevinSolver() { shutdown_thread_pool(); }

    // Article fragment: f_plus = 1, f_minus = alpha
    void set_two_level_f(double f_plus, double f_minus) {
        f_plus_ = f_plus;
        f_minus_ = f_minus;
    }

    LangevinEnsembleResult solve_ensemble(
        std::vector<DichotomicNoise>& noises,
        std::size_t N,
        std::size_t n_particles,
        double x0 = 0.0,
        std::size_t burn_in = 0,
        bool store_trajectory = true) {

        if (noises.size() < n_particles) {
            throw std::invalid_argument("noises.size() must be >= n_particles");
        }

        LangevinEnsembleResult res;
        res.L = L_;
        res.n_particles = n_particles;
        res.has_trajectory = store_trajectory;

        if (store_trajectory) {
            res.t.reserve(N);
            res.mean_x.reserve(N);
            res.mean_x_wrapped.reserve(N);
        }

        // wrapped координаты (для силы/потенциала), всегда в [0, L)
        std::vector<double> x(n_particles);
        // unwrapped координаты (для корректного <x(t)> и сравнения с теорией)
        std::vector<double> xu(n_particles);

        const double x0w = wrap_pos_(x0);
        for (std::size_t p = 0; p < n_particles; ++p) {
            x[p] = x0w;
            xu[p] = x0w;
        }

        // генераторы шума по потокам
        std::vector<std::mt19937> thread_rngs(n_threads_);
        for (std::size_t t = 0; t < n_threads_; ++t) {
            thread_rngs[t].seed(rng_() + static_cast<unsigned int>(t));
        }
        std::vector<std::normal_distribution<double>> thread_gaussians(
            n_threads_, std::normal_distribution<double>(0.0, 1.0));

        // статистика (по UNWRAPPED mean)
        bool has_prev_mean = false;
        double prev_mean_xu = 0.0;

        bool burnin_captured = false;
        double mean_xu_at_burnin = 0.0;

        double sum_msd = 0.0;
        std::size_t cnt_msd = 0;

        const double inv_zeta = D_ * beta_; // 1/ζ

        // Основной цикл
        for (std::size_t n = 0; n < N; ++n) {
            const std::size_t chunk = (n_particles + n_threads_ - 1) / n_threads_;

            std::vector<double> partial_sums_wrapped(n_threads_, 0.0);
            std::vector<double> partial_sums_unwrapped(n_threads_, 0.0);

            std::vector<std::atomic<int>> done(n_threads_);
            for (auto& d : done) d.store(0, std::memory_order_relaxed);

            std::size_t n_active = 0;

            for (std::size_t t = 0; t < n_threads_; ++t) {
                const std::size_t p_begin = t * chunk;
                const std::size_t p_end = std::min(n_particles, p_begin + chunk);
                if (p_begin >= p_end) break;

                ++n_active;

                workers_[t]->submit_task([this,
                                         &x,
                                         &xu,
                                         &noises,
                                         &partial_sums_wrapped,
                                         &partial_sums_unwrapped,
                                         &thread_rngs,
                                         &thread_gaussians,
                                         &done,
                                         inv_zeta,
                                         t,
                                         p_begin,
                                         p_end]() {
                    double local_sum_w = 0.0;
                    double local_sum_u = 0.0;

                    auto& local_rng = thread_rngs[t];
                    auto& gauss = thread_gaussians[t];

                    for (std::size_t p = p_begin; p < p_end; ++p) {
                        // дихотомика: f(t) через знак sigma
                        const double sigma = noises[p].step();
                        const double f_t = compute_f_(sigma); // f_+=1, f_- = alpha

                        // белый шум w_n ~ N(0,1)
                        const double w_n = gauss(local_rng);

                        // ---- предиктор x^* ----
                        const double dUdx_n = f_t * dVdx_(x[p]);   // U'(x,t)=f(t)V'(x)
                        const double drift_n = -inv_zeta * dUdx_n; // -ζ^{-1} U'
                        const double noise_term = sqrt2D_coeff_ * w_n * sqrt_dt_;
                        const double x_tilde = x[p] + drift_n * dt_ + noise_term;

                        // ---- корректор ----
                        const double dUdx_tilde = f_t * dVdx_(x_tilde);
                        const double drift_corr = -inv_zeta * 0.5 * (dUdx_n + dUdx_tilde);

                        const double delta = drift_corr * dt_ + noise_term;

                        // Обновляем UNWRAPPED координату
                        xu[p] += delta;

                        // А wrapped координату получаем из unwrapped
                        x[p] = wrap_pos_(xu[p]);

                        local_sum_w += x[p];
                        local_sum_u += xu[p];
                    }

                    partial_sums_wrapped[t] = local_sum_w;
                    partial_sums_unwrapped[t] = local_sum_u;
                    done[t].store(1, std::memory_order_release);
                });
            }

            // ждём потоки
            for (std::size_t t = 0; t < n_active; ++t) {
                while (done[t].load(std::memory_order_acquire) == 0) {
                    std::this_thread::yield();
                }
            }

            // среднее по ансамблю
            double sum_xw = 0.0;
            double sum_xu = 0.0;
            for (std::size_t t = 0; t < n_active; ++t) {
                sum_xw += partial_sums_wrapped[t];
                sum_xu += partial_sums_unwrapped[t];
            }

            const double mean_x_wrapped_n = sum_xw / static_cast<double>(n_particles);
            const double mean_x_unwrapped_n = sum_xu / static_cast<double>(n_particles);

            if (store_trajectory) {
                res.t.push_back(static_cast<double>(n) * dt_);
                res.mean_x.push_back(mean_x_unwrapped_n);
                res.mean_x_wrapped.push_back(mean_x_wrapped_n);
            }

            // статистика по mean_x_unwrapped_n
            if (!has_prev_mean) {
                prev_mean_xu = mean_x_unwrapped_n;
                has_prev_mean = true;
            } else {
                const double dmean = mean_x_unwrapped_n - prev_mean_xu;
                prev_mean_xu = mean_x_unwrapped_n;

                if (n > burn_in && n > 0) {
                    sum_msd += dmean * dmean;
                    ++cnt_msd;
                }
            }

            if (!burnin_captured && n >= burn_in) {
                mean_xu_at_burnin = mean_x_unwrapped_n;
                burnin_captured = true;
            }

            // скорость (по unwrapped mean)
            const double t_max = (N > burn_in) ? (static_cast<double>(N - burn_in) * dt_) : 0.0;
            if (burnin_captured && t_max > 0.0) {
                res.mean_velocity = (mean_x_unwrapped_n - mean_xu_at_burnin) / t_max;
            } else {
                res.mean_velocity = 0.0;
            }
        }

        // эффективный коэффициент диффузии (как и раньше: по инкрементам mean)
        if (cnt_msd > 0) {
            const double mean_msd = sum_msd / static_cast<double>(cnt_msd);
            res.diffusion_coeff = mean_msd / (2.0 * dt_);
        } else {
            res.diffusion_coeff = 0.0;
        }

        return res;
    }

private:
    void validate_params() const {
        if (dt_ <= 0.0) throw std::invalid_argument("dt must be > 0");
        if (L_ <= 0.0) throw std::invalid_argument("L must be > 0");
        if (D_ <= 0.0) throw std::invalid_argument("D must be > 0");
        if (beta_ <= 0.0) throw std::invalid_argument("beta must be > 0");
    }

    void initialize_thread_pool() {
        workers_.clear();
        workers_.reserve(n_threads_);
        for (std::size_t i = 0; i < n_threads_; ++i) {
            workers_.push_back(std::make_unique<WorkerThread>());
        }
    }

    void shutdown_thread_pool() {
        for (auto& w : workers_) {
            if (w) w->stop();
        }
    }

    inline double compute_f_(double sigma) const {
        if (mod_type_ == ModulationType::CONSTANT) return 1.0;
        return (sigma >= 0.0) ? f_plus_ : f_minus_;
    }

    inline double wrap_pos_(double v) const {
        v = std::fmod(v, L_);
        if (v < 0.0) v += L_;
        return v;
    }

    DVDXFunc dVdx_;
    double L_ = 1.0;
    double dt_ = 1e-3;

    ModulationType mod_type_ = ModulationType::CONSTANT;

    double D_ = 1.0;
    double beta_ = 1.0;

    double f_plus_ = 1.0;
    double f_minus_ = 0.0;

    std::mt19937 rng_;

    double sqrt2D_coeff_ = 1.0;
    double sqrt_dt_ = 0.0;
    double L_half_ = 0.0;

    std::size_t n_threads_ = 4;
    std::vector<std::unique_ptr<WorkerThread>> workers_;
};

#endif
