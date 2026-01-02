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
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"

// ============================================================================
// Структура для одной траектории
// ============================================================================
struct LangevinTrajectory {
    std::vector<double> t;
    std::vector<double> x;
    double L;
    double mean_position;
    double mean_velocity;
    double diffusion_coeff;
};

// ============================================================================
// Структура для ансамбля частиц
// ============================================================================
struct LangevinEnsembleResult {
    std::vector<double> t;
    std::vector<double> mean_x;
    double L;
    double mean_velocity;
    double diffusion_coeff;
    std::size_t n_particles;
    bool has_trajectory;
};

// ============================================================================
// Рабочий поток в пуле
// ============================================================================
class WorkerThread {
public:
    WorkerThread() : should_stop_(false), is_working_(false), is_running_(true) {
        thread_ = std::thread(&WorkerThread::run, this);
    }

    // Запрет копирования
    WorkerThread(const WorkerThread&) = delete;
    WorkerThread& operator=(const WorkerThread&) = delete;

    // Разрешение перемещения
    WorkerThread(WorkerThread&& other) noexcept
        : thread_(std::move(other.thread_)),
          task_queue_(std::move(other.task_queue_)),
          should_stop_(other.should_stop_),
          is_working_(other.is_working_),
          is_running_(other.is_running_) {
        // Сбрасываем другой объект в безопасное состояние
        other.should_stop_ = true;
        other.is_working_ = false;
        other.is_running_ = false;
    }

    WorkerThread& operator=(WorkerThread&& other) noexcept {
        if (this != &other) {
            stop(); // Останавливаем свой поток перед перемещением

            thread_ = std::move(other.thread_);
            task_queue_ = std::move(other.task_queue_);
            should_stop_ = other.should_stop_;
            is_working_ = other.is_working_;
            is_running_ = other.is_running_;

            other.should_stop_ = true;
            other.is_working_ = false;
            other.is_running_ = false;
        }
        return *this;
    }

    ~WorkerThread() {
        stop();
    }

    void submit_task(std::function<void()> task) {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (!is_running_) return; // Не добавляем в остановленный поток
            task_queue_.push(std::move(task));
        }
        cv_.notify_one();
    }

    void stop() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (!is_running_) return; // Уже остановлен
            should_stop_ = true;
            is_running_ = false;
        }
        cv_.notify_one();
        if (thread_.joinable()) {
            thread_.join();
        }
    }

    void wait_idle() {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [this]() { return task_queue_.empty() && !is_working_; });
    }

private:
    void run() {
        while (true) {
            std::function<void()> task;

            {
                std::unique_lock<std::mutex> lock(mutex_);
                cv_.wait(lock, [this]() { return !task_queue_.empty() || should_stop_; });

                if (should_stop_ && task_queue_.empty()) {
                    break;
                }

                if (!task_queue_.empty()) {
                    task = std::move(task_queue_.front());
                    task_queue_.pop();
                    is_working_ = true;
                }
            }

            if (task) {
                task();

                {
                    std::lock_guard<std::mutex> lock(mutex_);
                    is_working_ = false;
                }
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
    bool is_running_; // Флаг что поток запущен и не был остановлен
};

// ============================================================================
// Решатель Ланжевена с пулом потоков
// ============================================================================
class LangevinSolver {
public:
    using DVDXFunc = std::function<double(double)>;

    enum class ModulationType {
        CONSTANT,
        DICHOTOM_SYMMETRIC,
        EPSILON_PLUS_DICHOTOM
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
          sqrt2_coeff_(std::sqrt(0.2)),
          L_half_(L / 2.0),
          n_threads_(std::thread::hardware_concurrency()) {

        if (n_threads_ == 0) n_threads_ = 4;
        if (dt_ <= 0.0) throw std::invalid_argument("dt must be > 0");
        if (L_ <= 0.0) throw std::invalid_argument("L must be > 0");
        if (V0_ < 0.0) throw std::invalid_argument("V0 must be >= 0");

        // Инициализируем пул рабочих потоков
        initialize_thread_pool();
    }

    ~LangevinSolver() {
        shutdown_thread_pool();
    }

    // ========================================================================
    // Старый метод: одна частица
    // ========================================================================
    LangevinTrajectory solve(DichotomicNoise& noise,
                             std::size_t N,
                             double x0 = 0.0,
                             double epsilon = 0.0,
                             std::size_t burn_in = 0) {
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

        double sum_x = 0.0;
        std::size_t cnt = 0;
        for (std::size_t n = burn_in; n < N; ++n) {
            sum_x += res.x[n];
            ++cnt;
        }

        res.mean_position = (cnt > 0) ? (sum_x / static_cast<double>(cnt)) : 0.0;

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
    // Ансамбль с общим шумом (оптимизированная версия с корректной синхронизацией)
    // ========================================================================
    LangevinEnsembleResult solve_ensemble(DichotomicNoise& noise,
                                          std::size_t N,
                                          std::size_t n_particles,
                                          double x0 = 0.0,
                                          double epsilon = 0.0,
                                          std::size_t burn_in = 0,
                                          bool store_trajectory = true) {
        LangevinEnsembleResult res;
        res.L = L_;
        res.n_particles = n_particles;
        res.has_trajectory = store_trajectory;

        if (store_trajectory) {
            res.t.reserve(N);
            res.mean_x.reserve(N);
        }

        std::vector<DichotomicNoise> noise_copies;
        for (std::size_t t = 0; t < n_threads_; ++t) {
            noise_copies.push_back(noise);
        }

        std::vector<double> x(n_particles);
        for (auto& xi : x) {
            xi = std::fmod(x0, L_);
            if (xi < 0.0) xi += L_;
        }

        std::vector<std::mt19937> thread_rngs(n_threads_);
        for (std::size_t t = 0; t < n_threads_; ++t) {
            thread_rngs[t].seed(rng_() + t);
        }

        std::vector<std::normal_distribution<double>> thread_gaussians(
            n_threads_, std::normal_distribution<double>(0.0, 1.0));

        double sum_msd = 0.0;
        std::size_t cnt_msd = 0;
        double prev_mean_x = 0.0;

        for (std::size_t n = 0; n < N; ++n) {
            double sigma_n = noise_copies[0].step();
            double f_t = compute_f_t(sigma_n, epsilon);

            std::size_t chunk = (n_particles + n_threads_ - 1) / n_threads_;
            std::vector<double> partial_sums(n_threads_, 0.0);
            std::vector<std::atomic<int>> tasks_done(n_threads_);
            for (auto& td : tasks_done) {
                td.store(0, std::memory_order_relaxed);
            }

            std::size_t n_active = 0;

            // Отправляем задачи на все потоки
            for (std::size_t t = 0; t < n_threads_; ++t) {
                std::size_t p_begin = t * chunk;
                std::size_t p_end = std::min(n_particles, p_begin + chunk);

                if (p_begin >= p_end) break;

                ++n_active;

                workers_[t].submit_task([this, &x, &partial_sums, &thread_rngs,
                                        &thread_gaussians, &tasks_done, f_t, t,
                                        p_begin, p_end]() {
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

                    partial_sums[t] = local_sum;
                    tasks_done[t].store(1, std::memory_order_release);
                });
            }

            // Ждём завершения только активных задач
            for (std::size_t t = 0; t < n_active; ++t) {
                std::size_t p_begin = t * chunk;
                std::size_t p_end = std::min(n_particles, p_begin + chunk);
                if (p_begin >= p_end) continue;

                // Спин-ожидание с yield
                while (tasks_done[t].load(std::memory_order_acquire) == 0) {
                    std::this_thread::yield();
                }
            }

            double sum_x = 0.0;
            for (std::size_t t = 0; t < n_active; ++t) {
                sum_x += partial_sums[t];
            }
            double mean_x_n = sum_x / static_cast<double>(n_particles);

            if (store_trajectory) {
                res.t.push_back(static_cast<double>(n) * dt_);
                res.mean_x.push_back(mean_x_n);
            }

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

        if (cnt_msd > 0) {
            double mean_msd = sum_msd / static_cast<double>(cnt_msd);
            res.diffusion_coeff = mean_msd / (2.0 * dt_);
        } else {
            res.diffusion_coeff = 0.0;
        }

        return res;
    }

    // ========================================================================
    // Ансамбль с независимыми шумами (оптимизированная версия)
    // ========================================================================
    LangevinEnsembleResult solve_ensemble_independent(
        std::vector<DichotomicNoise>& noises,
        std::size_t N,
        std::size_t n_particles,
        double x0 = 0.0,
        double epsilon = 0.0,
        std::size_t burn_in = 0,
        bool store_trajectory = true) {

        LangevinEnsembleResult res;
        res.L = L_;
        res.n_particles = n_particles;
        res.has_trajectory = store_trajectory;

        if (store_trajectory) {
            res.t.reserve(N);
            res.mean_x.reserve(N);
        }

        std::vector<double> x(n_particles);
        for (auto& xi : x) {
            xi = std::fmod(x0, L_);
            if (xi < 0.0) xi += L_;
        }

        std::vector<std::mt19937> thread_rngs(n_threads_);
        for (std::size_t t = 0; t < n_threads_; ++t) {
            thread_rngs[t].seed(rng_() + t);
        }

        std::vector<std::normal_distribution<double>> thread_gaussians(
            n_threads_, std::normal_distribution<double>(0.0, 1.0));

        double sum_msd = 0.0;
        std::size_t cnt_msd = 0;
        double prev_mean_x = 0.0;

        for (std::size_t n = 0; n < N; ++n) {
            std::size_t chunk = (n_particles + n_threads_ - 1) / n_threads_;
            std::vector<double> partial_sums(n_threads_, 0.0);
            std::vector<std::atomic<int>> tasks_done(n_threads_);
            for (auto& td : tasks_done) {
                td.store(0, std::memory_order_relaxed);
            }

            std::size_t n_active = 0;

            for (std::size_t t = 0; t < n_threads_; ++t) {
                std::size_t p_begin = t * chunk;
                std::size_t p_end = std::min(n_particles, p_begin + chunk);

                if (p_begin >= p_end) break;

                ++n_active;

                workers_[t].submit_task([this, &x, &noises, &partial_sums,
                                        &thread_rngs, &thread_gaussians,
                                        &tasks_done, epsilon, t, p_begin, p_end]() {
                    double local_sum = 0.0;
                    auto& local_rng = thread_rngs[t];
                    auto& local_gaussian = thread_gaussians[t];

                    for (std::size_t p = p_begin; p < p_end; ++p) {
                        double sigma_n = noises[p].step();
                        double f_t = compute_f_t(sigma_n, epsilon);
                        double zeta_n = local_gaussian(local_rng);
                        double dVdx_n = dVdx_(x[p]);

                        double dxdt = -V0_ * dVdx_n * f_t + sqrt2_coeff_ * zeta_n;
                        x[p] += dxdt * dt_;
                        x[p] = std::fmod(x[p], L_);
                        if (x[p] < 0.0) x[p] += L_;

                        local_sum += x[p];
                    }

                    partial_sums[t] = local_sum;
                    tasks_done[t].store(1, std::memory_order_release);
                });
            }

            // Ждём завершения только активных задач
            for (std::size_t t = 0; t < n_active; ++t) {
                std::size_t p_begin = t * chunk;
                std::size_t p_end = std::min(n_particles, p_begin + chunk);
                if (p_begin >= p_end) continue;

                while (tasks_done[t].load(std::memory_order_acquire) == 0) {
                    std::this_thread::yield();
                }
            }

            double sum_x = 0.0;
            for (std::size_t t = 0; t < n_active; ++t) {
                sum_x += partial_sums[t];
            }
            double mean_x_n = sum_x / static_cast<double>(n_particles);

            if (store_trajectory) {
                res.t.push_back(static_cast<double>(n) * dt_);
                res.mean_x.push_back(mean_x_n);
            }

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

        if (cnt_msd > 0) {
            double mean_msd = sum_msd / static_cast<double>(cnt_msd);
            res.diffusion_coeff = mean_msd / (2.0 * dt_);
        } else {
            res.diffusion_coeff = 0.0;
        }

        return res;
    }

    // Геттеры
    double L() const { return L_; }
    double dt() const { return dt_; }
    double V0() const { return V0_; }
    std::size_t num_threads() const { return n_threads_; }

private:
    DVDXFunc dVdx_;
    double V0_, L_, dt_;
    ModulationType mod_type_;
    double sqrt2_coeff_;
    double L_half_;
    std::mt19937 rng_;
    std::normal_distribution<double> gaussian_;
    std::size_t n_threads_;
    std::vector<WorkerThread> workers_;

    void initialize_thread_pool() {
        workers_.clear();
        workers_.reserve(n_threads_);
        for (std::size_t i = 0; i < n_threads_; ++i) {
            workers_.emplace_back();
        }
    }

    void shutdown_thread_pool() {
        for (auto& w : workers_) {
            w.stop();
        }
    }

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