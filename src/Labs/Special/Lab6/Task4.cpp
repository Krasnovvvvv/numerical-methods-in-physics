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
#include <algorithm>

using namespace fdtd;

struct PhCResult {
    std::vector<double> lam_nm;
    std::vector<double> R, T;
    int periods;
    std::string label;
};

PhCResult runPhotonicCrystal(int num_periods,
                             double L_SiO2_nm,
                             double L_TiO2_nm,
                             const std::string& label,
                             int num_steps = 120000) {
    std::cout << "  PhC: " << num_periods << " periods, "
              << "SiO2=" << L_SiO2_nm << " nm, TiO2=" << L_TiO2_nm << " nm\n";

    const double c0       = 3e8;
    const double lam_min  = 300e-9;
    const double lam_max  = 900e-9;
    const double n_SiO2   = 1.45;
    const double n_TiO2   = 2.28;
    const double eps_SiO2 = n_SiO2 * n_SiO2;
    const double eps_TiO2 = n_TiO2 * n_TiO2;

    const double L_SiO2 = L_SiO2_nm * 1e-9;
    const double L_TiO2 = L_TiO2_nm * 1e-9;
    const double period = L_SiO2 + L_TiO2;

    const double dx   = 5e-9;
    const double Q    = 0.5;
    const int pml_N   = 40;

    const int left_gap_cells  = 140;
    const int right_gap_cells = 140;
    const int src_pos         = pml_N + 35;

    const double struct_len   = num_periods * period;
    const int struct_cells    = static_cast<int>(std::ceil(struct_len / dx));
    const int struct_start    = pml_N + left_gap_cells;
    const int struct_end      = struct_start + struct_cells;

    const int mon_ref   = src_pos + (struct_start - src_pos) / 2;
    const int mon_trans = struct_end + 60;
    const int N_total   = mon_trans + right_gap_cells + pml_N;

    auto freqs = freqsForRange(lam_min, lam_max, 700, c0);

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

    double x_pos = struct_start * dx;
    for (int p = 0; p < num_periods; ++p) {
        sim.addLayer({x_pos, x_pos + L_SiO2, eps_SiO2, 1.0, 0.0});
        x_pos += L_SiO2;
        sim.addLayer({x_pos, x_pos + L_TiO2, eps_TiO2, 1.0, 0.0});
        x_pos += L_TiO2;
    }

    sim.setSource(src_pos,
                  makeGaussianPulseForRange(lam_min, lam_max, c0, 6.0),
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
    res.R.resize(freqs.size());
    res.T.resize(freqs.size());

    for (size_t k = 0; k < freqs.size(); ++k) {
        res.lam_nm[k] = rt.wavelengths[k] * 1e9;
        res.R[k] = std::max(0.0, rt.R[k]);
        res.T[k] = std::max(0.0, rt.T[k]);
    }

    return res;
}

int main() {
    std::cout << "=== Task 4: Photonic Crystal ===\n\n";

    const double n_SiO2 = 1.45;
    const double n_TiO2 = 2.28;
    const double L_SiO2_opt_nm = 100.0;
    const double L_TiO2_opt_nm = L_SiO2_opt_nm * n_SiO2 / n_TiO2; // ~63.6 nm

    std::vector<int> periods_list = {5, 10, 20};

    {
        std::cout << "[1] Equal physical thickness (100/100 nm)\n";
        std::vector<PhCResult> results;
        for (int np : periods_list) {
            results.push_back(runPhotonicCrystal(np, 100.0, 100.0,
                std::to_string(np) + " periods"));
        }

        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);

        for (auto& r : results) {
            auto p = plot(r.lam_nm, r.R);
            p->display_name("R, " + r.label);
            p->line_width(1.8);
        }

        xlabel("Wavelength (nm)");
        ylabel("R");
        matplot::title("PhC (SiO2 100 nm / TiO2 100 nm) — Reflection");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.05});
        show();

        auto fig2 = figure(true);
        fig2->size(1000, 600);
        hold(on);

        auto& r20 = results.back();
        auto p1 = plot(r20.lam_nm, r20.R);
        p1->display_name("R");
        p1->line_width(2.0);

        auto p2 = plot(r20.lam_nm, r20.T);
        p2->display_name("T");
        p2->line_width(2.0);

        std::vector<double> RT_sum(r20.R.size());
        for (size_t k = 0; k < RT_sum.size(); ++k) {
            RT_sum[k] = r20.R[k] + r20.T[k];
        }

        auto p3 = plot(r20.lam_nm, RT_sum);
        p3->display_name("R+T");
        p3->line_style("--");
        p3->line_width(1.2);
        p3->color("gray");

        xlabel("Wavelength (nm)");
        ylabel("Coefficient");
        matplot::title("PhC 100/100 nm, 20 periods — R, T, R+T");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.05});
        show();
    }

    std::vector<PhCResult> results_optical;
    {
        std::cout << "\n[2] Equal optical thickness (100 / "
                  << L_TiO2_opt_nm << " nm)\n";

        for (int np : periods_list) {
            results_optical.push_back(runPhotonicCrystal(
                np, L_SiO2_opt_nm, L_TiO2_opt_nm,
                std::to_string(np) + " periods"));
        }

        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);

        for (auto& r : results_optical) {
            auto p = plot(r.lam_nm, r.R);
            p->display_name("R, " + r.label);
            p->line_width(1.8);
        }

        xlabel("Wavelength (nm)");
        ylabel("R");
        matplot::title("PhC (equal optical thickness) — Reflection");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.05});
        show();

        auto fig2 = figure(true);
        fig2->size(1000, 600);
        hold(on);

        auto& r20 = results_optical.back();
        auto p1 = plot(r20.lam_nm, r20.R);
        p1->display_name("R");
        p1->line_width(2.0);

        auto p2 = plot(r20.lam_nm, r20.T);
        p2->display_name("T");
        p2->line_width(2.0);

        std::vector<double> RT_sum(r20.R.size());
        for (size_t k = 0; k < RT_sum.size(); ++k) {
            RT_sum[k] = r20.R[k] + r20.T[k];
        }

        auto p3 = plot(r20.lam_nm, RT_sum);
        p3->display_name("R+T");
        p3->line_style("--");
        p3->line_width(1.2);
        p3->color("gray");

        xlabel("Wavelength (nm)");
        ylabel("Coefficient");
        matplot::title("Equal optical thickness, 20 periods — R, T, R+T");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.05});
        show();
    }

    {
        std::cout << "\n[3] Comparison: equal physical vs equal optical (20 periods)\n";

        auto r_phys = runPhotonicCrystal(20, 100.0, 100.0, "100/100");

        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);

        auto p1 = plot(r_phys.lam_nm, r_phys.R);
        p1->display_name("R: 100/100 nm (equal physical)");
        p1->line_width(2.0);

        auto p2 = plot(results_optical.back().lam_nm, results_optical.back().R);
        p2->display_name("R: equal optical thickness");
        p2->line_width(2.0);

        xlabel("Wavelength (nm)");
        ylabel("R");
        matplot::title("Comparison: equal physical vs equal optical thickness, 20 periods");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.05});
        show();
    }

    std::cout << "\nTask 4 complete.\n";
}