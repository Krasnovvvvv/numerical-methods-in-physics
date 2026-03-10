#ifndef NUMERICAL_METHODS_IN_PHYSICS_DIMLESSLANGEVINSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_DIMLESSLANGEVINSOLVER_H

#pragma once

#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <functional>
#include <memory>
#include <mutex>
#include <queue>
#include <random>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"
#include "Labs/Special/Additional/Base/RatchetParams.h"

struct EnsembleResult {
    std::vector<double> sim_time;
    std::size_t n_particles = 0;
    bool has_trajectory = false;

    std::vector<double> mean_x;
    double mean_velocity = 0.0;
};

class WorkerThread {
    void run() {
        while (true) {
            std::function<void()> task;

            {
                std::unique_lock<std::mutex> lock(mutex_);
                cv_.wait(lock, [this]() {
                    return !task_queue_.empty() || should_stop_;
                });
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
    bool is_running_;
    bool is_working_;

public:
    WorkerThread()
        : should_stop_(false), is_running_(true), is_working_(false) {
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
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (!is_running_) return;
            should_stop_ = true;
            is_running_ = false;
        }
        cv_.notify_one();
        if (thread_.joinable()) thread_.join();
    }
};

class DimLessLangevinSolver {
    void initialize_thread_pool() {
        workers_.clear();
        workers_.reserve(n_threads_);
        for (std::size_t i = 0; i < n_threads_; ++i) {
            workers_.push_back(std::make_unique<WorkerThread>());
        }
    }

    void shutdown_thread_pool() const {
        for (auto& w : workers_) {
            if (w) w->stop();
        }
    }

    [[nodiscard]] inline double compute_f(double sigma) const {
        return (sigma >= 0.0) ? 1.0 : params_.alpha;
    }

    RatchetParams params_;
    std::mt19937 rng_;
    std::size_t n_threads_ = 2;
    std::vector<std::unique_ptr<WorkerThread>> workers_;

public:
    explicit DimLessLangevinSolver(
        const RatchetParams &params,
        unsigned int seed = std::random_device{}())
        : params_(params),
          rng_(seed),
          n_threads_(std::thread::hardware_concurrency()) {
        if (n_threads_ == 0) n_threads_ = 4;
        initialize_thread_pool();
    }

    ~DimLessLangevinSolver() { shutdown_thread_pool(); }

EnsembleResult solve(
    std::vector<DichotomicNoise>& noises,
    std::size_t n_particles,
    std::size_t total_time,
    std::size_t burn_in = 0,
    double x0 = 0.0,
    bool store_trajectory = true) {

    if (noises.size() < n_particles) {
        throw std::invalid_argument("noises.size() must be >= n_particles");
    }

    EnsembleResult result;
    result.n_particles = n_particles;
    result.has_trajectory = store_trajectory;
    result.mean_velocity = 0.0;

    std::vector<double> x(n_particles, x0);

    const std::size_t N = static_cast<std::size_t>(total_time / params_.dt);
    if (burn_in == 0) burn_in = N / 10;

    if (store_trajectory) {
        result.sim_time.reserve(N);
        result.mean_x.reserve(N);
    }

    std::vector<std::mt19937> thread_rngs(n_threads_);
    for (std::size_t t = 0; t < n_threads_; ++t) {
        thread_rngs[t].seed(rng_() + static_cast<unsigned int>(t));
    }

    std::vector<std::normal_distribution<double>> thread_gaussians(
        n_threads_, std::normal_distribution<double>(0.0, 1.0));

    for (std::size_t n = 0; n < N; ++n) {
        const std::size_t chunk = (n_particles + n_threads_ - 1) / n_threads_;
        std::vector<double> partial_sums(n_threads_, 0.0);
        std::vector<std::atomic<int>> done(n_threads_);
        for (auto& d : done) d.store(0, std::memory_order_relaxed);

        std::size_t n_active = 0;

        for (std::size_t t = 0; t < n_threads_; ++t) {
            const std::size_t p_begin = t * chunk;
            const std::size_t p_end = std::min(n_particles, p_begin + chunk);
            if (p_begin >= p_end) break;

            ++n_active;

            workers_[t]->submit_task(
                [this, &x, &noises, &partial_sums, &thread_rngs, &thread_gaussians,
                 &done, t, p_begin, p_end]() {
                    double local_sum = 0.0;
                    auto& local_rng = thread_rngs[t];
                    auto& gauss = thread_gaussians[t];

                    for (std::size_t p = p_begin; p < p_end; ++p) {
                        const double sigma = noises[p].step();
                        const double f_t = compute_f(sigma);
                        const double w_n = gauss(local_rng);
                        const double noise = std::sqrt(2.0 * params_.dt) * w_n;

                        const double dUdx_n = f_t * params_.dVdx(x[p]);
                        const double x_tilde = x[p] - dUdx_n * params_.dt + noise;
                        const double dUdx_tilde = f_t * params_.dVdx(x_tilde);
                        const double delta =
                            -0.5 * (dUdx_n + dUdx_tilde) * params_.dt + noise;

                        x[p] += delta;
                        local_sum += x[p];
                    }

                    partial_sums[t] = local_sum;
                    done[t].store(1, std::memory_order_release);
                });
        }

        for (std::size_t t = 0; t < n_active; ++t) {
            while (done[t].load(std::memory_order_acquire) == 0) {
                std::this_thread::yield();
            }
        }

        double sum_x = 0.0;
        for (std::size_t t = 0; t < n_active; ++t) {
            sum_x += partial_sums[t];
        }

        const double mean_x = sum_x / static_cast<double>(n_particles);

        if (store_trajectory) {
            result.sim_time.emplace_back(static_cast<double>(n) * params_.dt);
            result.mean_x.emplace_back(mean_x);
        }
    }

    if (!store_trajectory || result.sim_time.size() < 2 || result.mean_x.size() < 2) {
        return result;
    }

    const std::size_t i0 = std::min(burn_in, result.sim_time.size() - 1);
    const std::size_t n_fit = result.sim_time.size() - i0;

    if (n_fit < 2) {
        return result;
    }

    double sum_t = 0.0;
    double sum_xm = 0.0;

    for (std::size_t i = i0; i < result.sim_time.size(); ++i) {
        sum_t += result.sim_time[i];
        sum_xm += result.mean_x[i];
    }

    const double mean_t = sum_t / static_cast<double>(n_fit);
    const double mean_xm = sum_xm / static_cast<double>(n_fit);

    double s_tt = 0.0;
    double s_tx = 0.0;

    for (std::size_t i = i0; i < result.sim_time.size(); ++i) {
        const double dt = result.sim_time[i] - mean_t;
        const double dx = result.mean_x[i] - mean_xm;
        s_tt += dt * dt;
        s_tx += dt * dx;
    }

    result.mean_velocity = (s_tt > 0.0) ? (s_tx / s_tt) : 0.0;
    return result;
}
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_DIMLESSLANGEVINSOLVER_H
