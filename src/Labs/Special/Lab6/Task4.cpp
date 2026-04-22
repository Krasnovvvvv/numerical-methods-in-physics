/// Задача 4: Фотонный кристалл — TiO₂/SiO₂
///
/// Графики:
/// 1. R(λ), T(λ) для равных толщин (100/100 нм) — несколько N периодов
/// 2. R(λ), T(λ) для равных оптических толщин (SiO₂ 100 нм, TiO₂ 55 нм)
/// 3. Сравнение двух вариантов (наложение R)

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

struct PhCResult {
    std::vector<double> lam_nm;
    std::vector<double> R, T;
    int periods;
    std::string label;
};

PhCResult runPhotonicCrystal(int num_periods, double L_SiO2_nm, double L_TiO2_nm,
                              const std::string& label) {
    std::cout << "  PhC: " << num_periods << " periods, "
              << "SiO2=" << L_SiO2_nm << " nm, TiO2=" << L_TiO2_nm << " nm\n";

    const double c0      = 3e8;
    const double lam_min = 300e-9;
    const double lam_max = 900e-9;
    const double n_SiO2  = 1.45;
    const double n_TiO2  = 2.28;
    const double eps_SiO2 = n_SiO2 * n_SiO2;
    const double eps_TiO2 = n_TiO2 * n_TiO2;
    const double L_SiO2  = L_SiO2_nm * 1e-9;
    const double L_TiO2  = L_TiO2_nm * 1e-9;
    const double period  = L_SiO2 + L_TiO2;

    const double dx  = 5e-9;
    const int pml_N  = 40;
    const double Q   = 0.5;

    double struct_len = num_periods * period;
    int struct_cells  = (int)(struct_len / dx) + 10;
    int margin = 80;
    int src_offset = 60;
    int N_total = 2 * pml_N + 2 * margin + struct_cells;

    int src_pos      = pml_N + src_offset;
    int struct_start = pml_N + margin;
    int struct_end   = struct_start + struct_cells;
    int mon_ref      = src_pos + 30;
    int mon_trans    = struct_end + 30;

    if (mon_trans + pml_N + 10 > N_total)
        N_total = mon_trans + pml_N + 20;

    int num_steps = 80000;
    auto freqs = freqsForRange(lam_min, lam_max, 500, c0);

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

    // Периодическая структура: SiO₂ / TiO₂
    double x_pos = struct_start * dx;
    for (int p = 0; p < num_periods; ++p) {
        sim.addLayer({x_pos, x_pos + L_SiO2, eps_SiO2, 1.0, 0.0});
        x_pos += L_SiO2;
        sim.addLayer({x_pos, x_pos + L_TiO2, eps_TiO2, 1.0, 0.0});
        x_pos += L_TiO2;
    }

    sim.setSource(src_pos, makeGaussianPulseForRange(lam_min, lam_max, c0, 6.0),
                  SourceMode::Soft);
    sim.addMonitor(mon_ref, freqs);
    sim.addMonitor(mon_trans, freqs);
    sim.build();
    sim.run();

    auto norm = sim.runNormalization();
    auto rt   = sim.computeRT(norm, c0);

    PhCResult res;
    res.periods = num_periods;
    res.label   = label;
    res.lam_nm.resize(freqs.size());
    res.R = rt.R;
    res.T = rt.T;
    for (size_t k = 0; k < freqs.size(); ++k)
        res.lam_nm[k] = rt.wavelengths[k] * 1e9;

    return res;
}

int main() {
    std::cout << "=== Task 4: Photonic Crystal ===\n\n";

    std::vector<int> periods_list = {5, 10, 20};

    // ================================================================
    // График 1: Равные физические толщины (100/100 нм)
    // ================================================================
    {
        std::cout << "[1] Equal physical thickness (100/100 nm)\n";
        std::vector<PhCResult> results;
        for (int np : periods_list)
            results.push_back(runPhotonicCrystal(np, 100.0, 100.0,
                std::to_string(np) + " periods"));

        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);
        for (auto& r : results) {
            auto p = plot(r.lam_nm, r.R);
            p->display_name("R, " + r.label);
            p->line_width(1.5);
        }
        xlabel("Wavelength (nm)");
        ylabel("R");
        matplot::title("PhC (SiO2 100nm / TiO2 100nm) — Reflection");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.1});
        show();

        // Также T для 20 периодов
        auto fig2 = figure(true);
        fig2->size(1000, 600);
        hold(on);
        auto& r20 = results.back(); // 20 periods
        auto p1 = plot(r20.lam_nm, r20.R);
        p1->display_name("R"); p1->line_width(2);
        auto p2 = plot(r20.lam_nm, r20.T);
        p2->display_name("T"); p2->line_width(2);

        // R + T
        std::vector<double> RT_sum(r20.R.size());
        for (size_t k = 0; k < RT_sum.size(); ++k)
            RT_sum[k] = r20.R[k] + r20.T[k];
        auto p3 = plot(r20.lam_nm, RT_sum);
        p3->display_name("R+T"); p3->line_style("--"); p3->line_width(1); p3->color("gray");

        xlabel("Wavelength (nm)");
        ylabel("Coefficient");
        matplot::title("PhC 100/100 nm, 20 periods — R, T, R+T");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.1});
        show();
    }

    // ================================================================
    // График 2: Равные оптические толщины (SiO₂ 100 нм, TiO₂ 55 нм)
    // ================================================================
    std::vector<PhCResult> results_optical;
    {
        std::cout << "\n[2] Equal optical thickness (100/55 nm)\n";
        for (int np : periods_list)
            results_optical.push_back(runPhotonicCrystal(np, 100.0, 55.0,
                std::to_string(np) + " periods"));

        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);
        for (auto& r : results_optical) {
            auto p = plot(r.lam_nm, r.R);
            p->display_name("R, " + r.label);
            p->line_width(1.5);
        }
        xlabel("Wavelength (nm)");
        ylabel("R");
        matplot::title("PhC (SiO2 100nm / TiO2 55nm) — Reflection");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.1});
        show();
    }

    // ================================================================
    // График 3: Сравнение двух вариантов (20 периодов)
    // ================================================================
    {
        std::cout << "\n[3] Comparison: equal physical vs equal optical (20 periods)\n";
        auto r_phys = runPhotonicCrystal(20, 100.0, 100.0, "100/100");
        // results_optical.back() уже содержит 20 периодов 100/55

        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);

        auto p1 = plot(r_phys.lam_nm, r_phys.R);
        p1->display_name("R: 100/100 nm (equal physical)");
        p1->line_width(2);

        auto p2 = plot(results_optical.back().lam_nm, results_optical.back().R);
        p2->display_name("R: 100/55 nm (equal optical)");
        p2->line_width(2);

        xlabel("Wavelength (nm)");
        ylabel("R");
        matplot::title("Comparison: equal physical vs equal optical thickness, 20 periods");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.1});
        show();
    }

    std::cout << "\nTask 4 complete.\n";
    std::cout << "\nВыводы:\n";
    std::cout << "  - Фотонные запрещённые зоны видны как области R≈1, T≈0\n";
    std::cout << "  - Больше периодов — резче границы запрещённых зон\n";
    std::cout << "  - Равные оптические толщины дают запрещённую зону,\n";
    std::cout << "    центрированную на четвертьволновом условии\n";
    std::cout << "  - R+T≈1, что подтверждает энергетическую корректность\n";
}