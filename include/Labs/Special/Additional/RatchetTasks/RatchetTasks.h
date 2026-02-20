#ifndef NUMERICAL_METHODS_IN_PHYSICS_RATCHETTASKS_H
#define NUMERICAL_METHODS_IN_PHYSICS_RATCHETTASKS_H

#pragma once

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "Labs/Special/Additional/DffusionSolvers/LangevinSolver.h"
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"
#include "Helpers/Plotter.h"
#include "Helpers/Timer.h"

class RatchetTasks {

    static constexpr double TWO_PI = 2.0 * M_PI;
    static constexpr double FOUR_PI = 4.0 * M_PI;

    // V'(x) для симметричного потенциала
    static auto getDVdxSymmetric() {
        return [](double x) -> double {
            return -std::cos(TWO_PI * x);
        };
    }

    // V'(x) для асимметричного потенциала
    static auto getDVdxAsymmetric() {
        return [](double x) -> double {
            return (std::cos(TWO_PI * x) + 0.5 * std::cos(FOUR_PI * x));
        };
    }

    static std::size_t choose_burn_in(std::size_t N) {
        return N / 10;
    }

public:

    /*

    static void taskA(double L, double dt, std::size_t N) {

        double a = 1.0;
        double tau_c = 0.05;
        std::size_t n_particles = 20000;
        std::vector<double> V0_vals = {0.8, 1.5, 2.5};
        std::vector<std::vector<double>> t_vecs;
        std::vector<std::vector<double>> x_vecs;
        std::vector<std::string> labels;

        std::size_t burn_in = choose_burn_in(N);
        double t_burn = burn_in * dt;
        double T_window = 100.0;
        double T_plot_start = t_burn;
        double T_plot_end = t_burn + T_window;

        for (double V0 : V0_vals) {
            LangevinSolver solver(
                getDVdxSymmetric(), V0, L, dt,
                LangevinSolver::ModulationType::DICHOTOM_SYMMETRIC, 42);

            // Создаём вектор независимых дихотомных шумов
            std::vector<DichotomicNoise> noises;
            noises.reserve(n_particles);
            for (std::size_t p = 0; p < n_particles; ++p) {
                noises.emplace_back(a, tau_c, dt,
                    42u + static_cast<unsigned int>(p));
            }

            auto res = solver.solve_ensemble_independent(
                noises, N, n_particles, 0.0,
                0.0, burn_in, true);

            std::vector<double> t_cut, x_cut;
            for (std::size_t i = 0; i < res.mean_x.size(); ++i) {
                double t = res.t[i];
                if (t < T_plot_start) continue;
                if (t > T_plot_end) break;
                t_cut.push_back(t);
                x_cut.push_back(res.mean_x[i]);
            }

            if (!t_cut.empty()) {
                t_vecs.push_back(std::move(t_cut));
                x_vecs.push_back(std::move(x_cut));
                labels.push_back("V₀ = " + std::to_string(V0));
            }
        }

        Plotter plotter;
        plotter.plot(t_vecs, x_vecs, labels, "t", "<x>");
    }

    static void taskB(double L, double dt, std::size_t N) {

        double a = 1.0;
        double tau_c = 1.9;
        std::size_t n_particles = 100000;
        std::vector<double> V0_vals = {0.8, 1.5, 2.5, 3.5};
        std::vector<std::vector<double>> t_vecs;
        std::vector<std::vector<double>> x_vecs;
        std::vector<std::string> labels;

        std::size_t burn_in = choose_burn_in(N);
        double t_burn = burn_in * dt;
        double T_window = 100.0;
        double T_plot_start = t_burn;
        double T_plot_end = t_burn + T_window;

        for (double V0 : V0_vals) {
            LangevinSolver solver(
                getDVdxAsymmetric(), V0, L, dt,
                LangevinSolver::ModulationType::DICHOTOM_SYMMETRIC, 200);

            // Вектор независимых шумов
            std::vector<DichotomicNoise> noises;
            noises.reserve(n_particles);
            for (std::size_t p = 0; p < n_particles; ++p) {
                noises.emplace_back(a, tau_c, dt, 200u + static_cast<unsigned int>(p));
            }

            auto res = solver.solve_ensemble_independent(noises, N, n_particles, 0.0, 0.0, burn_in, true);

            std::vector<double> t_cut, x_cut;
            for (std::size_t i = 0; i < res.mean_x.size(); ++i) {
                double t = res.t[i];
                if (t < T_plot_start) continue;
                if (t > T_plot_end) break;
                t_cut.push_back(t);
                x_cut.push_back(res.mean_x[i]);
            }

            if (!t_cut.empty()) {
                t_vecs.push_back(std::move(t_cut));
                x_vecs.push_back(std::move(x_cut));
                labels.push_back("V₀ = " + std::to_string(V0));
            }

            std::cout << V0 << "\t\t" << res.mean_velocity << "\n";
        }

        Plotter plotter;
        plotter.plot(t_vecs, x_vecs, labels, "t", "<x>");
    }

    static void taskC(double L, double dt, std::size_t N) {

        double a = 1.0;
        double V0 = 0.8;
        std::size_t n_particles = 100000;
        std::vector<double> tau_c_vals = {
            3.5, 3.0, 2.5, 2.0,
            1.9, 1.8, 1.7, 1.6,
            1.4, 1.0, 0.8};
        std::vector<double> inv_tau_c_vals;
        std::vector<double> v_mean_vals;

        std::size_t burn_in = choose_burn_in(N);

        for (double tau_c : tau_c_vals) {
            LangevinSolver solver(
                getDVdxAsymmetric(), V0, L, dt,
                LangevinSolver::ModulationType::DICHOTOM_SYMMETRIC, 300);

            // Вектор независимых шумов
            std::vector<DichotomicNoise> noises;
            noises.reserve(n_particles);
            for (std::size_t p = 0; p < n_particles; ++p) {
                noises.emplace_back(a, tau_c, dt,
                    300u + static_cast<unsigned int>(p));
            }

            auto res = solver.solve_ensemble_independent(noises,
                N, n_particles, 0.0, 0.0, burn_in, true);

            double inv_tau_c = 1.0 / tau_c;
            inv_tau_c_vals.push_back(inv_tau_c);
            v_mean_vals.push_back(res.mean_velocity);

        }

        Plotter plotter;
        plotter.plot(inv_tau_c_vals, v_mean_vals,
            "1/τ_c", "1/τ_c", "<v>");
    }

    static void taskG(double L, double dt, std::size_t N) {

        double gamma_b = 0.7;
        double V0 = 0.2;
        std::size_t n_particles = 100000;
        std::vector<double> epsilon_vals = {0.35};
        std::vector<double> v_mean_vals;
        std::vector<double> D_eff_vals;

        std::size_t burn_in = choose_burn_in(N);
        double t_burn = burn_in * dt;
        double T_window = 10.0;
        double T_plot_start = t_burn;
        double T_plot_end = t_burn + T_window;


        std::vector<double> t_demo, x_demo;

        for (double epsilon : epsilon_vals) {
            LangevinSolver solver(
                getDVdxAsymmetric(), V0, L, dt,
                LangevinSolver::ModulationType::EPSILON_PLUS_DICHOTOM, 400);

            // Вектор независимых шумов для каждой частицы
            std::vector<DichotomicNoise> noises;
            noises.reserve(n_particles);
            for (std::size_t p = 0; p < n_particles; ++p) {
                noises.push_back(
                    DichotomicNoise::ZeroMeanAsymmetric(2.0, -1.0, gamma_b, dt,
                                                        400u + static_cast<unsigned int>(p)));
            }

            bool store_traj = true;
            auto res = solver.solve_ensemble_independent(noises, N, n_particles, 0.0, epsilon, burn_in, store_traj);

            v_mean_vals.push_back(res.mean_velocity);
            D_eff_vals.push_back(res.diffusion_coeff);

            // Сохраняем первое окно для демонстрации
            if (t_demo.empty() && res.has_trajectory) {
                for (std::size_t i = 0; i < res.mean_x.size(); ++i) {
                    double t = res.t[i];
                    if (t < T_plot_start) continue;
                    if (t > T_plot_end) break;
                    t_demo.push_back(t);
                    x_demo.push_back(res.mean_x[i]);
                }
            }
        }

        Plotter plotter;

        if (!t_demo.empty()) {
            plotter.plot(t_demo, x_demo, "x(t)", "t", "");
        }

        plotter.plot(epsilon_vals, v_mean_vals, "(ε)", "ε", "");

    }

    static void taskD(double L, double dt, std::size_t N) {

    double gamma_b = 0.5;
    double V0 = 0.5;              // Амплитуда потенциала
    double u = 0.7;               // Средний уровень модуляции
    std::size_t n_particles = 100000;

    // Сетка по амплитуде w (амплитуда пульсации)
    std::vector<double> w_vals = {0.4};

    std::vector<double> v_mean_vals;
    std::vector<double> D_eff_vals;

    std::size_t burn_in = choose_burn_in(N);
    double t_burn = burn_in * dt;
    double T_window = 100.0;
    double T_plot_start = t_burn;
    double T_plot_end = t_burn + T_window;


    std::vector<double> t_demo, x_demo;

    for (double w : w_vals) {

        LangevinSolver solver(
            getDVdxAsymmetric(), V0, L, dt,
            LangevinSolver::ModulationType::AMPLITUDE_MODULATION, 500);

        // Вектор независимых симметричных дихотомных шумов
        // σ(t) = ±1 с одинаковыми вероятностями
        std::vector<DichotomicNoise> noises;
        noises.reserve(n_particles);

        for (std::size_t p = 0; p < n_particles; ++p) {
            noises.emplace_back(
                3.0,
                -5.0,
                1.9,
                0.6,
                dt,
                500u + static_cast<unsigned int>(p)
            );
        }

        bool store_traj = true;

        auto res = solver.solve_ensemble_amplitude(
            noises, N, n_particles, 0.0, u, w, burn_in, store_traj);

        v_mean_vals.push_back(res.mean_velocity);
        D_eff_vals.push_back(res.diffusion_coeff);

        std::cout << w << "\t\t\t" << res.mean_velocity
                  << "\t\t" << res.diffusion_coeff << "\n";

        // Сохраняем первую траекторию для демонстрации
        if (t_demo.empty() && res.has_trajectory) {
            for (std::size_t i = 0; i < res.mean_x.size(); ++i) {
                double t = res.t[i];
                if (t < T_plot_start) continue;
                if (t > T_plot_end) break;
                t_demo.push_back(t);
                x_demo.push_back(res.mean_x[i]);
            }
        }
    }

    Plotter plotter;
    if (!t_demo.empty()) {
        plotter.plot(t_demo, x_demo, "x(t)", "t", "");
    }

    plotter.plot(w_vals, v_mean_vals,
        "<v>(w)", "w (амплитуда)", "<v>");

}
    */

static void taskComparisonAnalyticalPaper(double dt, std::size_t N) {
    const double L = 2.0;

    const double V1 = 3.5;
    const double V2 = 0.5 * V1;

    const double alpha  = -1./3.0;
    const double f_plus = 1.0;
    const double f_minus = alpha;

    // термопараметры
    const double T  = 15.0;
    const double kB = 1.0;
    const double beta = 1.0 / (kB * T);

    // свободная диффузия
    const double D = 2.0; // kBT/friction

    // дихотомика
    const double a = 1.0;   // sigma = ±1
    const double tau_c = 0.08;

    const double tau_plus  = tau_c;
    const double tau_minus = tau_c;
    const double delta = tau_minus / (tau_plus + tau_minus); // 1/2

    const auto n_particles = static_cast<std::size_t>(1e4);
    const std::size_t burn_in = 20'000;

    // dV/dx для V(x)=V1 sin(2πx/L)+V2 sin(4πx/L)
    auto dVdx = [=](double x) {
        const double k1 = 2.0 * M_PI / L;
        const double k2 = 4.0 * M_PI / L;
        return V1 * k1 * std::cos(k1 * x) + V2 * k2 * std::cos(k2 * x);
    };

    const double two_pi = 2.0 * M_PI;
    const double xi = (L / two_pi) * (L / two_pi) / (D * 2.0 * tau_c);

    const double z = xi / (4.0 * delta * (1.0 - delta));
    const double denom = (1.0 + 4.0 * z) * (1.0 + 4.0 * z) * (1.0 + z);

    const double Phi1_S = (3.0 * xi * (1.0 + 2.0 * z)) / denom;
    const double Phi2_S = (1.0 - 2.0 * delta) * (6.0 * xi * z) / denom;

    const double prefactor = (M_PI * D) / (4.0 * L);
    const double v_analytical =
        prefactor * (beta * beta * beta) * (V1 * V1) * V2
        * (1.0 - alpha) * (1.0 - alpha)
        * ((1.0 + alpha) * Phi1_S + (1.0 - alpha) * Phi2_S);

    std::cout << std::setprecision(12)
              << "Parameters:\n"
              << "  L = " << L << "\n"
              << "  dt = " << dt << "\n"
              << "  N = " << N << "\n"
              << "  T_total = " << (static_cast<double>(N) * dt) << "\n"
              << "  T = " << T << "  beta=" << beta << "\n"
              << "  D (free) = " << D << "\n"
              << "  V1 = " << V1 << "  V2 = " << V2 << "\n"
              << "  alpha (f_-)= " << alpha << " (f_+=1)\n"
              << "аналитика скорости: " << v_analytical << "\n";

    // ---- численно ----
    LangevinSolver solver(dVdx, L, dt, D, beta, LangevinSolver::ModulationType::TWO_LEVEL_F, 600);
    solver.set_two_level_f(f_plus, f_minus);

    std::vector<DichotomicNoise> noises;
    noises.reserve(n_particles);
    for (std::size_t p = 0; p < n_particles; ++p) {
        noises.emplace_back(a, tau_c, dt, 600u + static_cast<unsigned>(p));
    }

    Timer<std::chrono::duration<double>> timer;
    auto res = solver.solve_ensemble(noises, N, n_particles, 0.0, burn_in, true);
    const double elapsed_sec = timer.elapsed();

    // mean_x теперь уже UNWRAPPED (важно для сравнения с v*t)
    const auto& mean_x_unwrapped = res.mean_x;

    // full-window v по UNWRAPPED mean_x
    double v_full = 0.0;
    if (res.t.size() >= 2 && mean_x_unwrapped.size() == res.t.size()) {
        const double dt_full = res.t.back() - res.t.front();
        if (dt_full > 0.0) {
            v_full = (mean_x_unwrapped.back() - mean_x_unwrapped.front()) / dt_full;
        }
    }

    const double error_abs = std::abs(v_full - v_analytical);
    const double error_rel = (v_analytical != 0.0)
        ? (error_abs / std::abs(v_analytical)) * 100.0
        : 100.0;

    std::cout << std::setprecision(12)
              << "v_numerical (solver.mean_velocity) = " << res.mean_velocity << "\n"
              << "v_numerical (full-window) = " << v_full << "\n"
              << "abs error = " << error_abs << "\n"
              << "rel error = " << error_rel << " %\n"
              << "simulation time = " << elapsed_sec << " sec (" << elapsed_sec / 60.0 << " min)\n";

    // ---- графики ----
    // x(t): численно vs линия аналитики (только после burn-in)
    const std::size_t step = std::max<std::size_t>(1, N / 800);

    std::vector<double> t_plot, x_plot, x_an;
    t_plot.reserve(N / step + 1);
    x_plot.reserve(N / step + 1);
    x_an.reserve(N / step + 1);

    const std::size_t i0 = std::min<std::size_t>(
        burn_in,
        (res.t.size() > 0 ? (res.t.size() - 1) : 0)
    );

    const double t0  = res.t.empty() ? 0.0 : res.t[i0];
    const double x0u = mean_x_unwrapped.empty() ? 0.0 : mean_x_unwrapped[i0];

    for (std::size_t i = i0; i < res.t.size() && i < mean_x_unwrapped.size(); i += step) {
        const double t = res.t[i];
        t_plot.push_back(t);
        x_plot.push_back(mean_x_unwrapped[i]);            // <-- корректные unwrapped координаты
        x_an.push_back(x0u + v_analytical * (t - t0));    // <-- теоретическая прямая
    }

    {
        std::vector<std::vector<double>> xs = {t_plot, t_plot};
        std::vector<std::vector<double>> ys = {x_plot, x_an};
        std::vector<std::string> labels = {
            "Численное решение (ансамбль, unwrapped)",
            "Аналитическое решение (формула 12)"
        };
        Plotter plotter;
        plotter.plot(xs, ys, labels, "Время t (сек)", "Позиция x(t)");
    }

    // v(t) = (x(t)-x(t0))/(t-t0), тоже считаем только после burn-in
    std::vector<double> t_v, v_est, v_an;
    for (std::size_t k = 5; k < t_plot.size(); ++k) {
        const double tt = t_plot[k];
        const double dt_rel = tt - t0;
        if (dt_rel <= 0.0) continue;

        t_v.push_back(tt);
        v_est.push_back((x_plot[k] - x0u) / dt_rel);
        v_an.push_back(v_analytical);
    }

    {
        std::vector<std::vector<double>> xs = {t_v, t_v};
        std::vector<std::vector<double>> ys = {v_est, v_an};
        std::vector<std::string> labels = {
            "(<x>-<x>0)/(t-t0)",
            "v_analytical"
        };
        Plotter plotter;
        plotter.plot(xs, ys, labels, "Время t (сек)", "Скорость v(t)");
    }
}

};

#endif // NUMERICAL_METHODS_IN_PHYSICS_RATCHETTASKS_H