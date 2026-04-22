/// Задача 3: Отражение от диэлектрической пластины (SiO₂)
///
/// Графики:
/// 1. R(λ), T(λ) при L=200 нм: FDTD vs аналитика Фабри-Перо
/// 2. R(λ) при разных L (100, 200, 400 нм)
/// 3. R(λ) при разном времени моделирования (5000, 15000, 30000, 60000 шагов)

#include "Labs/Special/Lab6/Base/Grid.h"
#include "Labs/Special/Lab6/Base/Source.h"
#include "Labs/Special/Lab6/Solver/FDTDSolver.h"
#include "Labs/Special/Lab6/Base/PML.h"
#include "Labs/Special/Lab6/Base/Monitor.h"
#include "Labs/Special/Lab6/Base/Material.h"
#include "Labs/Special/Lab6/Simulation/Simulation.h"

#include <matplot/matplot.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

using namespace fdtd;

/// Аналитика Фабри-Перо для пластины с n, толщиной L, при длине волны lambda
struct FPResult { double R, T; };

FPResult fabryPerot(double n, double L, double lambda) {
    double r12 = (1.0 - n) / (1.0 + n);
    double delta = 2.0 * M_PI * n * L / lambda;
    double r12_sq = r12 * r12;
    double denom = 1.0 + r12_sq * r12_sq - 2.0 * r12_sq * std::cos(2.0 * delta);
    double R = 2.0 * r12_sq * (1.0 - std::cos(2.0 * delta)) / denom;
    double T = std::pow(1.0 - r12_sq, 2) / denom;
    return {R, T};
}

/// Результат одного прогона
struct SlabResult {
    std::vector<double> lam_nm;
    std::vector<double> R, T;
    std::vector<double> R_th, T_th;
    double L_nm;
    int steps;
};

SlabResult runSlab(double L_nm, int num_steps) {
    const double c0        = 3e8;
    const double lambda_min = 300e-9;
    const double lambda_max = 900e-9;
    const double n_quartz   = 1.45;
    const double eps_quartz = n_quartz * n_quartz;
    const double L          = L_nm * 1e-9;

    const double dx = 5e-9;
    const int pml_N = 40;
    const int N_domain = 600;
    const int N_total  = N_domain + 2 * pml_N;
    const double Q = 0.5;

    int src_pos      = pml_N + 60;
    int mon_ref      = pml_N + 120;
    int center       = pml_N + N_domain / 2;
    int slab_start_i = center - (int)(L / (2.0 * dx));
    int slab_end_i   = center + (int)(L / (2.0 * dx));
    int mon_trans    = slab_end_i + 60;

    auto freqs = freqsForRange(lambda_min, lambda_max, 500, c0);

    Simulation::Config cfg;
    cfg.num_cells     = N_total;
    cfg.dx            = dx;
    cfg.courant       = Q;
    cfg.c_speed       = c0;
    cfg.num_steps     = num_steps;
    cfg.pml_thickness = pml_N;
    cfg.pml_order     = 3;
    cfg.pml_delta     = 1e-8;

    Simulation sim(cfg);
    sim.addLayer({slab_start_i * dx, slab_end_i * dx, eps_quartz, 1.0, 0.0});
    sim.setSource(src_pos, makeGaussianPulseForRange(lambda_min, lambda_max, c0, 6.0),
                  SourceMode::Soft);
    sim.addMonitor(mon_ref, freqs);
    sim.addMonitor(mon_trans, freqs);
    sim.build();
    sim.run();

    auto norm = sim.runNormalization();
    auto rt   = sim.computeRT(norm, c0);

    // Аналитика
    SlabResult res;
    res.L_nm  = L_nm;
    res.steps = num_steps;
    res.lam_nm.resize(freqs.size());
    res.R = rt.R;
    res.T = rt.T;
    res.R_th.resize(freqs.size());
    res.T_th.resize(freqs.size());

    for (size_t k = 0; k < freqs.size(); ++k) {
        res.lam_nm[k] = rt.wavelengths[k] * 1e9;
        auto fp = fabryPerot(n_quartz, L, rt.wavelengths[k]);
        res.R_th[k] = fp.R;
        res.T_th[k] = fp.T;
    }
    return res;
}

int main() {
    std::cout << "=== Task 3: Dielectric Slab ===\n\n";

    // ================================================================
    // График 1: R(λ), T(λ) для L=200 нм — FDTD vs Fabry-Perot
    // ================================================================
    {
        std::cout << "[1] L=200 nm, 30000 steps: FDTD vs Fabry-Perot\n";
        auto r = runSlab(200.0, 30000);

        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);

        auto p1 = plot(r.lam_nm, r.R);
        p1->display_name("FDTD R"); p1->line_width(2);
        auto p2 = plot(r.lam_nm, r.R_th);
        p2->display_name("Theory R"); p2->line_style("--"); p2->line_width(2);
        auto p3 = plot(r.lam_nm, r.T);
        p3->display_name("FDTD T"); p3->line_width(2);
        auto p4 = plot(r.lam_nm, r.T_th);
        p4->display_name("Theory T"); p4->line_style("--"); p4->line_width(2);

        xlabel("Wavelength (nm)");
        ylabel("Coefficient");
        matplot::title("SiO2 slab L=200 nm — FDTD vs Fabry-Perot");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.1});
        show();
    }

    // ================================================================
    // График 2: R(λ) при разных L
    // ================================================================
    {
        std::cout << "\n[2] Varying slab thickness\n";
        std::vector<double> thicknesses = {100.0, 200.0, 400.0};

        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);

        for (double L : thicknesses) {
            std::cout << "  L=" << L << " nm...\n";
            auto r = runSlab(L, 30000);
            auto p = plot(r.lam_nm, r.R);
            p->display_name("L=" + std::to_string((int)L) + " nm");
            p->line_width(1.5);
        }

        xlabel("Wavelength (nm)");
        ylabel("R");
        matplot::title("Reflection from SiO2 slab — varying thickness");
        matplot::legend();
        xlim({300.0, 900.0});
        show();
    }

    // ================================================================
    // График 3: Влияние времени моделирования на R(λ)
    // ================================================================
    {
        std::cout << "\n[3] Varying simulation time\n";
        std::vector<int> step_counts = {50, 500, 1000};

        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);

        for (int steps : step_counts) {
            std::cout << "  steps=" << steps << "...\n";
            auto r = runSlab(200.0, steps);
            auto p = plot(r.lam_nm, r.R);
            p->display_name(std::to_string(steps) + " steps");
            p->line_width(1.5);
        }

        // Добавим теорию для сравнения
        auto r_ref = runSlab(200.0, 60000);
        auto p = plot(r_ref.lam_nm, r_ref.R_th);
        p->display_name("Theory (Fabry-Perot)");
        p->line_style("--"); p->line_width(2); p->color("black");

        xlabel("Wavelength (nm)");
        ylabel("R");
        matplot::title("Convergence with simulation time (L=200 nm)");
        matplot::legend();
        xlim({300.0, 900.0});
        show();
    }

    std::cout << "\nTask 3 complete.\n";
    std::cout << "\nВыводы:\n";
    std::cout << "  - R(λ) и T(λ) демонстрируют осцилляции Фабри-Перо\n";
    std::cout << "  - При увеличении L число осцилляций растёт\n";
    std::cout << "  - Малое время моделирования недостаточно: учитываются не все\n";
    std::cout << "    многократные отражения внутри пластины → спектр неточен\n";
    std::cout << "  - При достаточном числе шагов FDTD сходится к аналитике\n";
}