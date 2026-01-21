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

class RatchetTasks {
private:
    // Предвычисленные константы для потенциалов
    static constexpr double TWO_PI = 2.0 * M_PI;
    static constexpr double FOUR_PI = 4.0 * M_PI;

    // V'(x) для симметричного потенциала (задание A)
    static auto getDVdxSymmetric() {
        return [](double x) -> double {
            return -std::cos(TWO_PI * x);
        };
    }

    // V'(x) для асимметричного потенциала (задания B, C, G)
    static auto getDVdxAsymmetric() {
        return [](double x) -> double {
            return (std::cos(TWO_PI * x) + 0.5 * std::cos(FOUR_PI * x));
        };
    }

    static std::size_t choose_burn_in(std::size_t N) {
        return N / 10;
    }

public:
    // =========================
    // ЗАДАНИЕ А
    // =========================
    static void taskA(double L, double dt, std::size_t N) {
        std::cout << "\n>>> ЗАДАНИЕ А: x(t) во времени (симметричный потенциал)\n";
        std::cout << " V'(x) = cos(2πx)\n";
        std::cout << " f(t) = (1/2)(1 + σ(t))\n";
        std::cout << " ОЖИДАНИЕ: НЕТ линейного роста x(t)\n";
        std::cout << "───────────────────────────────────────────────────────────\n";

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
                noises.emplace_back(a, tau_c, dt, 42u + static_cast<unsigned int>(p));
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

            std::cout << "V₀ = " << V0
                      << ": средняя скорость ≈ " << res.mean_velocity
                      << " (ожидается ≈ 0)\n";
        }

        Plotter plotter;
        plotter.plot(t_vecs, x_vecs, labels, "t", "<x>");
        std::cout << "✓ График x(t) (A) отрисован!\n";
    }

    // =========================
    // ЗАДАНИЕ Б
    // =========================
    static void taskB(double L, double dt, std::size_t N) {
        std::cout << "\n>>> ЗАДАНИЕ Б: x(t) во времени (асимметричный потенциал)\n";
        std::cout << " V'(x) = cos(2πx) + (1/2)cos(4πx)\n";
        std::cout << " f(t) = (1/2)(1 + σ(t))\n";
        std::cout << " ОЖИДАНИЕ: x(t) с ненулевым наклоном\n";
        std::cout << "───────────────────────────────────────────────────────────\n";

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

        std::cout << "V₀\t\t\n";
        std::cout << "─────────────────────────────────────\n";

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
        std::cout << "✓ График x(t) (B) отрисован!\n";
    }

    // =========================
    // ЗАДАНИЕ В: v(1/τ_c)
    // =========================
    static void taskC(double L, double dt, std::size_t N) {
        std::cout << "\n>>> ЗАДАНИЕ В: v vs 1/τ_c (резонанс)\n";
        std::cout << " V'(x) = cos(2πx) + (1/2)cos(4πx)\n";
        std::cout << " f(t) = (1/2)(1 + σ(t))\n";
        std::cout << "───────────────────────────────────────────────────────────\n";

        double a = 1.0;
        double V0 = 0.8;
        std::size_t n_particles = 100000;
        std::vector<double> tau_c_vals = {3.5, 3.0, 2.5, 2.0, 1.9, 1.8, 1.7, 1.6, 1.4, 1.0, 0.8};
        std::vector<double> inv_tau_c_vals;
        std::vector<double> v_mean_vals;

        std::size_t burn_in = choose_burn_in(N);

        std::cout << "τ_c\t\t1/τ_c\t\t\n";
        std::cout << "──────────────────────────────────────────────\n";

        for (double tau_c : tau_c_vals) {
            LangevinSolver solver(
                getDVdxAsymmetric(), V0, L, dt,
                LangevinSolver::ModulationType::DICHOTOM_SYMMETRIC, 300);

            // Вектор независимых шумов
            std::vector<DichotomicNoise> noises;
            noises.reserve(n_particles);
            for (std::size_t p = 0; p < n_particles; ++p) {
                noises.emplace_back(a, tau_c, dt, 300u + static_cast<unsigned int>(p));
            }

            auto res = solver.solve_ensemble_independent(noises,
                N, n_particles, 0.0, 0.0, burn_in, true);

            double inv_tau_c = 1.0 / tau_c;
            inv_tau_c_vals.push_back(inv_tau_c);
            v_mean_vals.push_back(res.mean_velocity);

            std::cout << tau_c << "\t\t" << inv_tau_c
                      << "\t\t" << res.mean_velocity << "\n";
        }

        Plotter plotter;
        plotter.plot(inv_tau_c_vals, v_mean_vals, "1/τ_c", "1/τ_c", "<v>");
        std::cout << "✓ График v(1/τ_c) отрисован!\n";
    }

    // =========================
    // ЗАДАНИЕ Г: f(t) = ε + σ(t)
    // =========================
    static void taskG(double L, double dt, std::size_t N) {
        std::cout << "\n>>> ЗАДАНИЕ Г: Пульсирующий рэтчет f(t) = ε + σ(t)\n";
        std::cout << " Асимметричный дихотомный шум (±1, gamma_a != gamma_b)\n";
        std::cout << "───────────────────────────────────────────────────────────\n";

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

        std::cout << "ε\t\t\t\tD_eff\n";
        std::cout << "────────────────────────────────────────────────\n";

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

            std::cout << epsilon << "\t\t" << res.mean_velocity
                      << "\t\t" << res.diffusion_coeff << "\n";

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
            std::cout << "\n Отрисовка траектории (окно 100 c после burn-in)...\n";
            plotter.plot(t_demo, x_demo, "x(t)", "t", "");
        }

        std::cout << "\n Отрисовка графика (ε)...\n";
        plotter.plot(epsilon_vals, v_mean_vals, "(ε)", "ε", "");

        std::cout << " Отрисовка графика D_eff(ε)...\n";
        plotter.plot(epsilon_vals, D_eff_vals, "D_eff(ε)", "ε", "D_eff");

        std::cout << "✓ Графики для задания Г отрисованы!\n";
    }

    static void taskD(double L, double dt, std::size_t N) {

    std::cout << "\n>>> ЗАДАНИЕ D: Амплитудный рэтчет f(t) = u + w*σ(t)\n";
    std::cout << " Асимметричный потенциал с модуляцией АМПЛИТУДЫ\n";
    std::cout << " V'(x) = cos(2πx) + (1/2)cos(4πx)\n";
    std::cout << "───────────────────────────────────────────────────────────\n";

    double gamma_b = 0.5;          // → tau_c ≈ 0.67 с
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

    std::cout << "w (амплитуда)\t\t<v>\t\t\tD_eff\n";
    std::cout << "────────────────────────────────────────────────\n";

    std::vector<double> t_demo, x_demo;

    for (double w : w_vals) {

        LangevinSolver solver(
            getDVdxAsymmetric(), V0, L, dt,
            LangevinSolver::ModulationType::AMPLITUDE_MODULATION, 500);

        // Вектор независимых симметричных дихотомных шумов
        // σ(t) = ±1 с одинаковыми вероятностями
        std::vector<DichotomicNoise> noises;
        noises.reserve(n_particles);

        //double tau_c_noise = 1.0 / gamma_b;  // время корреляции шума

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

        // Вызов НОВОГО метода с амплитудной модуляцией
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
        std::cout << "\n Отрисовка траектории (w=0)...\n";
        plotter.plot(t_demo, x_demo, "x(t)", "t", "");
    }

    std::cout << "\n Отрисовка графика <v>(w)...\n";
    plotter.plot(w_vals, v_mean_vals, "<v>(w)", "w (амплитуда)", "<v>");

    std::cout << " Отрисовка графика D_eff(w)...\n";
    plotter.plot(w_vals, D_eff_vals, "D_eff(w)", "w (амплитуда)", "D_eff");

    std::cout << "✓ Графики для задания D (амплитудный рэтчет) отрисованы!\n";
}

// =========================
// СРАВНЕНИЕ С АНАЛИТИКОЙ (Формулы 12-13 из статьи)
// Rozenbaum et al. JETP Lett. 79, 388 (2004)
// =========================

static void taskComparisonAnalytical(double L, double dt, std::size_t N) {

    double L_period = 1.0;

    double V0 = 3.5;            // Масштаб потенциальной энергии
    double V1 = V0;             // Амплитуда первой гармоники
    double V2 = 0.5 * V0;       // Амплитуда второй гармоники

    // Температура
    double T = 15.0;            // Абсолютная температура (в условных единицах)
    double kB = 1.0;            // Постоянная Больцмана (естественные единицы)
    double kT = kB * T;         // Произведение
    double beta = 1.0 / kT;

    double a = 1.0;             // Амплитуда шума

    double tau_c = 1.0;         // τ_c — время корреляции

    double gamma_b = 0.8;  // Переходная вероятность
    double D = 1.0 / (beta * gamma_b);  // Коэффициент диффузии

    double tau_plus = tau_c;     // τ_+ = τ_c
    double tau_minus = tau_c;    // τ_- = τ_c
    double tau = 1.0 / (1.0/tau_plus + 1.0/tau_minus);  // Среднее время
    double delta = tau_minus / (tau_plus + tau_minus);

    std::size_t n_particles = 100000;
    std::size_t burn_in = 8'000;

    // ========================================================================
    // ВЫВОД ПАРАМЕТРОВ
    // ========================================================================

    std::cout << "\n┌─ ПАРАМЕТРЫ МОДЕЛИ ──────────────────────────────────────┐\n";
    std::cout << "│ ГЕОМЕТРИЯ ПОТЕНЦИАЛА:                                    │\n";
    std::cout << "│  L (период потенциала)               = " << std::setw(12) << L_period << "   │\n";
    std::cout << "│  V(x) = V₁·sin(2πx/L) + V₂·sin(4πx/L)                   │\n";
    std::cout << "│  V₀ (масштаб амплитуды)              = " << std::setw(12) << V0 << "   │\n";
    std::cout << "│  V₁ (первая гармоника)               = " << std::setw(12) << V1 << "   │\n";
    std::cout << "│  V₂ (вторая гармоника)               = " << std::setw(12) << V2 << "   │\n";

    std::cout << "│                                                          │\n";
    std::cout << "│ ТЕРМОДИНАМИКА:                                           │\n";
    std::cout << "│  T (абсолютная температура)          = " << std::setw(12) << T << "   │\n";
    std::cout << "│  k_B·T                               = " << std::setw(12) << kT << "   │\n";
    std::cout << "│  β = (k_B·T)⁻¹ (обратная темп.)      = " << std::setw(12) << beta << "   │\n";
    std::cout << "│  β·V₀ (безразм. амплитуда потенц.)   = " << std::setw(12) << (beta * V0) << "   │\n";
    std::cout << "│  β·V₀ << 1? (высокотемп. предел)     = ";
    std::cout << ((beta * V0 < 0.1) ? "ДА ✓" : "НЕТ ✗") << std::setw(24) << "│\n";

    std::cout << "│                                                          │\n";
    std::cout << "│ СТОХАСТИЧЕСКИЕ ПАРАМЕТРЫ:                               │\n";
    std::cout << "│  a (амплитуда шума σ = ±1)            = " << std::setw(12) << a << "   │\n";
    std::cout << "│  τ_c (время корреляции)               = " << std::setw(12) << tau_c << "   │\n";
    std::cout << "│  γ_- = 1/τ_c (вероятность)            = " << std::setw(12) << gamma_b << "   │\n";
    std::cout << "│  D = a²/(2γ_-) (коэфф. диффузии)      = " << std::setw(12) << D << "   │\n";
    std::cout << "│  τ_+ = τ_- = τ_c (симметричный)       = " << std::setw(12) << tau_plus << "   │\n";
    std::cout << "│  δ = τ_-(τ_+ + τ_-)⁻¹                 = " << std::setw(12) << delta << "   │\n";

    std::cout << "│                                                          │\n";
    std::cout << "│ ЧИСЛЕННОЕ МОДЕЛИРОВАНИЕ:                                │\n";
    std::cout << "│  Количество частиц в ансамбле         = " << std::setw(12) << n_particles << "   │\n";
    std::cout << "│  Общее число шагов                    = " << std::setw(12) << N << "   │\n";
    std::cout << "│  Шаг интегрирования dt                = " << std::setw(12) << dt << "   │\n";
    std::cout << "└──────────────────────────────────────────────────────────┘\n";

    double two_pi = 2.0 * M_PI;
    double xi = (L_period / two_pi) * (L_period / two_pi) / (D * tau_c);

    double z = xi / (4.0 * delta * (1.0 - delta));


    double denom = (1.0 + 4.0*z) * (1.0 + 4.0*z) * (1.0 + z);

    // Φ₁⁽ˢ⁾(ξ, δ) = 3ξ(1 + 2z) / [(1+4z)²(1+z)]
    double Phi1_S = (3.0 * xi * (1.0 + 2.0*z)) / denom;


    // Φ₂⁽ˢ⁾(ξ, δ) = (1 - 2δ) · 6ξz / [(1+4z)²(1+z)]
    double Phi2_S = (1.0 - 2.0*delta) * (6.0 * xi * z) / denom;

    std::cout << "\n┌─ АНАЛИТИЧЕСКОЕ РЕШЕНИЕ (Формула 12) ────────────────────┐\n";

    double alpha = (V2 - 0.0) / (V1 + V2);  // Асимметрия: (V₂)/(V₁ + V₂)

    // ν⁽ˢ'ᴰ⁾ = (π·D)/(4L) · β³ · V₁² · V₂ · (1-α)² · [(1+α)Φ₁⁽ˢ'ᴰ⁾ + (1-α)Φ₂⁽ˢ'ᴰ⁾]

    double prefactor = (M_PI * D) / (4.0 * L_period);

    double beta3 = beta * beta * beta;

    double V1_sq = V1 * V1;

    double asym_factor = (1.0 - alpha) * (1.0 - alpha);

    double Phi_weighted = (1.0 + alpha) * Phi1_S + (1.0 - alpha) * Phi2_S;

    double v_analytical = prefactor * beta3 * V1_sq * V2 * asym_factor * Phi_weighted;
            std::cout << "│ ν_АНАЛИТИЧЕСКАЯ = " << std::scientific << std::setprecision(8)
              << v_analytical << " │\n" << std::defaultfloat;


    LangevinSolver solver(
        getDVdxAsymmetric(), V0, L, dt,
        LangevinSolver::ModulationType::DICHOTOM_SYMMETRIC, 600);

    std::vector<DichotomicNoise> noises;
    noises.reserve(n_particles);
    for (std::size_t p = 0; p < n_particles; ++p) {
        noises.emplace_back(a, tau_c, dt, 600u + static_cast<unsigned>(p));
    }

    auto res_numerical = solver.solve_ensemble_independent(
        noises, N, n_particles, 0.0, 0.0, burn_in, true);

    std::cout << "│ Time window (burn-in period): " << std::setw(19) << (burn_in * dt) << " с  │\n";
    std::cout << "│ Time window (measurement):    " << std::setw(19)
              << ((N - burn_in) * dt) << " с  │\n";
    std::cout << "│                                                          │\n";
    std::cout << "│ v_ЧИСЛЕННАЯ = " << std::scientific << std::setprecision(8)
              << res_numerical.mean_velocity << " │\n" << std::defaultfloat;

    // ========================================================================
    // СРАВНЕНИЕ И АНАЛИЗ ОШИБКИ
    // ========================================================================

    double error_abs = std::abs(res_numerical.mean_velocity - v_analytical);
    double error_rel = (v_analytical != 0.0)
        ? (error_abs / std::abs(v_analytical)) * 100.0
        : 100.0;

    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║               СРАВНЕНИЕ РЕЗУЛЬТАТОВ                         ║\n";
    std::cout << "╠════════════════════════════════════════════════════════════╣\n";
    std::cout << "║ v_аналитическая = " << std::scientific << std::setprecision(8)
              << v_analytical << std::defaultfloat << "                     ║\n";
    std::cout << "║ v_численная     = " << std::scientific << std::setprecision(8)
              << res_numerical.mean_velocity << std::defaultfloat << "                     ║\n";
    std::cout << "║ Абсолютная ошибка = " << std::scientific << std::setprecision(4)
              << error_abs << std::defaultfloat << "                          ║\n";
    std::cout << "║ Относительная ошибка = " << std::fixed << std::setprecision(2)
              << error_rel << " %" << std::setw(38) << "║\n";
    std::cout << "╠════════════════════════════════════════════════════════════╣\n";

    if (error_rel < 5.0) {
        std::cout << "║ ✓✓✓ ОТЛИЧНОЕ СОГЛАСИЕ (ошибка < 5%)!                      ║\n";
    } else if (error_rel < 15.0) {
        std::cout << "║ ✓✓  ХОРОШЕЕ СОГЛАСИЕ (ошибка 5-15%)                       ║\n";
    } else if (error_rel < 30.0) {
        std::cout << "║ ✓   ПРИЕМЛЕМОЕ СОГЛАСИЕ (ошибка 15-30%)                   ║\n";
    } else if (error_rel < 50.0) {
        std::cout << "║ ≈   УМЕРЕННОЕ РАСХОЖДЕНИЕ (ошибка 30-50%)                 ║\n";
    } else {
        std::cout << "║ ✗   ЗНАЧИТЕЛЬНОЕ РАСХОЖДЕНИЕ (ошибка > 50%)!              ║\n";
        std::cout << "║ → Проверьте параметры и формулы!                          ║\n";
    }
    std::cout << "╚════════════════════════════════════════════════════════════╝\n";

    std::size_t burn_in_idx = burn_in;
    std::size_t T_plot_end = std::min(static_cast<std::size_t>(N), burn_in + 10000);
    std::size_t step = std::max<std::size_t>(1, (T_plot_end - burn_in_idx) / 500);

    std::vector<double> t_plot, v_numerical_plot, v_analytical_plot;

    double x0_window = res_numerical.mean_x[burn_in_idx];

    for (std::size_t i = burn_in_idx; i < T_plot_end; i += step) {
        if (i < res_numerical.t.size()) {
            t_plot.push_back(res_numerical.t[i]);
            v_numerical_plot.push_back(res_numerical.mean_x[i]);

            double t_rel = res_numerical.t[i] - res_numerical.t[burn_in_idx];
            v_analytical_plot.push_back(x0_window + v_analytical * t_rel);
        }
    }

    std::vector<std::vector<double>> xs = {t_plot, t_plot};
    std::vector<std::vector<double>> ys = {v_numerical_plot, v_analytical_plot};
    std::vector<std::string> labels = {
        "Численное решение (ансамбль)",
        "Аналитическое решение (формула 12)"
    };

    Plotter plotter;
    plotter.plot(xs, ys, labels, "Время t (сек)", "Позиция x(t)");

}
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_RATCHETTASKS_H