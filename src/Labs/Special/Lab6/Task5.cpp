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

struct CavityResult {
    std::vector<double> lam_nm;
    std::vector<double> R, T, RTsum;
    std::string label;
};

CavityResult runBraggCavity(int n_periods, double cavity_nm, double lambda_center_nm,
                            const std::string& label) {
    std::cout << "  Bragg cavity: " << n_periods << " per/side, cavity="
              << cavity_nm << " nm, lambda0=" << lambda_center_nm << " nm\n";

    const double c0       = 3e8;
    const double lam_min  = 300e-9;
    const double lam_max  = 900e-9;
    const double n_SiO2   = 1.45;
    const double n_TiO2   = 2.28;
    const double eps_SiO2 = n_SiO2 * n_SiO2;
    const double eps_TiO2 = n_TiO2 * n_TiO2;

    const double lam_c    = lambda_center_nm * 1e-9;
    const double L_SiO2   = lam_c / (4.0 * n_SiO2);
    const double L_TiO2   = lam_c / (4.0 * n_TiO2);
    const double L_cavity = cavity_nm * 1e-9;
    const double period   = L_SiO2 + L_TiO2;

    const double dx  = 5e-9;
    const int pml_N  = 40;
    const double Q   = 0.5;

    const double struct_len = 2.0 * n_periods * period + L_cavity;
    const int struct_cells  = static_cast<int>(std::ceil(struct_len / dx)) + 8;

    const int left_gap_cells  = 120;
    const int right_gap_cells = 120;
    const int src_pos         = pml_N + 30;
    const int struct_start    = pml_N + left_gap_cells;
    const int struct_end      = struct_start + struct_cells;
    const int mon_ref         = struct_start - 30;
    const int mon_trans       = struct_end + 30;
    const int N_total         = struct_end + right_gap_cells + pml_N;

    const int num_steps = 160000;
    auto freqs = freqsForRange(lam_min, lam_max, 800, c0);

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

    // left mirror: TiO2 / SiO2
    double x_pos = struct_start * dx;
    for (int p = 0; p < n_periods; ++p) {
        sim.addLayer({x_pos, x_pos + L_TiO2, eps_TiO2, 1.0, 0.0});
        x_pos += L_TiO2;
        sim.addLayer({x_pos, x_pos + L_SiO2, eps_SiO2, 1.0, 0.0});
        x_pos += L_SiO2;
    }

    // air cavity
    x_pos += L_cavity;

    // right mirror: SiO2 / TiO2
    for (int p = 0; p < n_periods; ++p) {
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

    CavityResult res;
    res.label = label;
    res.lam_nm.resize(freqs.size());
    res.R = rt.R;
    res.T = rt.T;
    res.RTsum.resize(freqs.size());

    for (size_t k = 0; k < freqs.size(); ++k) {
        res.lam_nm[k] = rt.wavelengths[k] * 1e9;
        res.RTsum[k]  = rt.R[k] + rt.T[k];
    }

    return res;
}

CavityResult runPurePhC(int total_periods, double lambda_center_nm) {
    std::cout << "  Pure PhC: " << total_periods << " periods\n";

    const double c0       = 3e8;
    const double lam_min  = 300e-9;
    const double lam_max  = 900e-9;
    const double n_SiO2   = 1.45;
    const double n_TiO2   = 2.28;
    const double eps_SiO2 = n_SiO2 * n_SiO2;
    const double eps_TiO2 = n_TiO2 * n_TiO2;

    const double lam_c    = lambda_center_nm * 1e-9;
    const double L_SiO2   = lam_c / (4.0 * n_SiO2);
    const double L_TiO2   = lam_c / (4.0 * n_TiO2);
    const double period   = L_SiO2 + L_TiO2;

    const double dx  = 5e-9;
    const int pml_N  = 40;
    const double Q   = 0.5;

    const double struct_len = total_periods * period;
    const int struct_cells  = static_cast<int>(std::ceil(struct_len / dx)) + 8;

    const int left_gap_cells  = 120;
    const int right_gap_cells = 120;
    const int src_pos         = pml_N + 30;
    const int struct_start    = pml_N + left_gap_cells;
    const int struct_end      = struct_start + struct_cells;
    const int mon_ref         = struct_start - 30;
    const int mon_trans       = struct_end + 30;
    const int N_total         = struct_end + right_gap_cells + pml_N;

    auto freqs = freqsForRange(lam_min, lam_max, 800, c0);

    Simulation::Config cfg;
    cfg.num_cells     = N_total;
    cfg.dx            = dx;
    cfg.courant       = Q;
    cfg.c_speed       = c0;
    cfg.num_steps     = 120000;
    cfg.pml_thickness = pml_N;
    cfg.pml_order     = 3;
    cfg.pml_delta     = 1e-8;

    Simulation sim(cfg);

    double x_pos = struct_start * dx;
    for (int p = 0; p < total_periods; ++p) {
        sim.addLayer({x_pos, x_pos + L_TiO2, eps_TiO2, 1.0, 0.0});
        x_pos += L_TiO2;
        sim.addLayer({x_pos, x_pos + L_SiO2, eps_SiO2, 1.0, 0.0});
        x_pos += L_SiO2;
    }

    sim.setSource(src_pos, makeGaussianPulseForRange(lam_min, lam_max, c0, 6.0),
                  SourceMode::Soft);
    sim.addMonitor(mon_ref, freqs);
    sim.addMonitor(mon_trans, freqs);
    sim.build();
    sim.run();

    auto norm = sim.runNormalization();
    auto rt   = sim.computeRT(norm, c0);

    CavityResult res;
    res.label = "Pure PhC (" + std::to_string(total_periods) + " per)";
    res.lam_nm.resize(freqs.size());
    res.R = rt.R;
    res.T = rt.T;
    res.RTsum.resize(freqs.size());

    for (size_t k = 0; k < freqs.size(); ++k) {
        res.lam_nm[k] = rt.wavelengths[k] * 1e9;
        res.RTsum[k]  = rt.R[k] + rt.T[k];
    }

    return res;
}

void runFieldSnapshotsCW(int n_periods, double cavity_nm, double lambda_center_nm) {
    std::cout << "  Field snapshots (CW): " << n_periods << " per/side\n";

    const double c0       = 3e8;
    const double n_SiO2   = 1.45;
    const double n_TiO2   = 2.28;
    const double eps_SiO2 = n_SiO2 * n_SiO2;
    const double eps_TiO2 = n_TiO2 * n_TiO2;

    const double lam_c    = lambda_center_nm * 1e-9;
    const double f0       = c0 / lam_c;
    const double L_SiO2   = lam_c / (4.0 * n_SiO2);
    const double L_TiO2   = lam_c / (4.0 * n_TiO2);
    const double L_cavity = cavity_nm * 1e-9;
    const double period   = L_SiO2 + L_TiO2;

    const double dx  = 5e-9;
    const int pml_N  = 40;
    const double Q   = 0.5;

    const double struct_len = 2.0 * n_periods * period + L_cavity;
    const int struct_cells  = static_cast<int>(std::ceil(struct_len / dx)) + 8;

    const int left_gap_cells  = 120;
    const int right_gap_cells = 120;
    const int src_pos         = pml_N + 30;
    const int struct_start    = pml_N + left_gap_cells;
    const int struct_end      = struct_start + struct_cells;
    const int N_total         = struct_end + right_gap_cells + pml_N;

    GridConfig gc;
    gc.num_cells = N_total;
    gc.dx        = dx;
    gc.courant   = Q;
    gc.c_speed   = c0;
    Grid grid(gc);

    std::vector<Layer> layers;
    double x_pos = struct_start * dx;

    for (int p = 0; p < n_periods; ++p) {
        layers.push_back({x_pos, x_pos + L_TiO2, eps_TiO2, 1.0, 0.0});
        x_pos += L_TiO2;
        layers.push_back({x_pos, x_pos + L_SiO2, eps_SiO2, 1.0, 0.0});
        x_pos += L_SiO2;
    }

    double cavity_x0 = x_pos;
    x_pos += L_cavity;
    double cavity_x1 = x_pos;

    for (int p = 0; p < n_periods; ++p) {
        layers.push_back({x_pos, x_pos + L_SiO2, eps_SiO2, 1.0, 0.0});
        x_pos += L_SiO2;
        layers.push_back({x_pos, x_pos + L_TiO2, eps_TiO2, 1.0, 0.0});
        x_pos += L_TiO2;
    }

    applyLayers(grid, layers);

    PMLConfig pml_cfg;
    pml_cfg.thickness = pml_N;
    pml_cfg.order     = 3;
    pml_cfg.delta     = 1e-8;
    applyPML(grid, pml_cfg);

    auto wf = makeCW(f0, 20.0 / f0, 4.0);
    Source source{src_pos, wf, SourceMode::Soft};

    const int num_steps = 120000;
    std::vector<int> snap_times = {80000, 90000, 100000, 110000};
    std::vector<std::vector<double>> snap_Ey;
    std::vector<int> snap_labels;

    std::vector<double> x_nm(N_total);
    for (int i = 0; i < N_total; ++i) x_nm[i] = grid.xE(i) * 1e9;

    int snap_idx = 0;
    for (int n = 0; n < num_steps; ++n) {
        FDTDSolver::updateE(grid);
        source.inject(grid, n);
        FDTDSolver::updateH(grid);

        if (snap_idx < (int)snap_times.size() && n == snap_times[snap_idx]) {
            std::vector<double> ey(N_total);
            for (int i = 0; i < N_total; ++i) ey[i] = grid.Ey[i];
            snap_Ey.push_back(ey);
            snap_labels.push_back(n);
            ++snap_idx;
        }
    }

    double max_abs = 1e-12;
    for (const auto& s : snap_Ey) {
        for (double v : s) max_abs = std::max(max_abs, std::abs(v));
    }

    using namespace matplot;
    auto fig = figure(true);
    fig->size(1200, 600);
    hold(on);

    for (size_t s = 0; s < snap_Ey.size(); ++s) {
        auto p = plot(x_nm, snap_Ey[s]);
        p->display_name("n=" + std::to_string(snap_labels[s]));
        p->line_width(1.3);
    }

    auto l1 = plot({cavity_x0 * 1e9, cavity_x0 * 1e9}, {-max_abs, max_abs});
    l1->line_style("--"); l1->color("gray"); l1->line_width(1.0);
    l1->display_name("Cavity edges");

    auto l2 = plot({cavity_x1 * 1e9, cavity_x1 * 1e9}, {-max_abs, max_abs});
    l2->line_style("--"); l2->color("gray"); l2->line_width(1.0);

    xlabel("x (nm)");
    ylabel("E_y");
    matplot::title("CW field localization in Bragg microcavity");
    matplot::legend();
    show();
}

void runFieldSnapshotsGaussian(int n_periods, double cavity_nm, double lambda_center_nm) {
    std::cout << "  Field snapshots (Gaussian): " << n_periods << " per/side\n";

    const double c0       = 3e8;
    const double lam_min  = 300e-9;
    const double lam_max  = 900e-9;
    const double n_SiO2   = 1.45;
    const double n_TiO2   = 2.28;
    const double eps_SiO2 = n_SiO2 * n_SiO2;
    const double eps_TiO2 = n_TiO2 * n_TiO2;

    const double lam_c    = lambda_center_nm * 1e-9;
    const double L_SiO2   = lam_c / (4.0 * n_SiO2);
    const double L_TiO2   = lam_c / (4.0 * n_TiO2);
    const double L_cavity = cavity_nm * 1e-9;
    const double period   = L_SiO2 + L_TiO2;

    const double dx  = 5e-9;
    const int pml_N  = 40;
    const double Q   = 0.5;

    const double struct_len = 2.0 * n_periods * period + L_cavity;
    const int struct_cells  = static_cast<int>(std::ceil(struct_len / dx)) + 8;

    const int left_gap_cells  = 120;
    const int right_gap_cells = 120;
    const int src_pos         = pml_N + 30;
    const int struct_start    = pml_N + left_gap_cells;
    const int struct_end      = struct_start + struct_cells;
    const int N_total         = struct_end + right_gap_cells + pml_N;

    const int num_steps = 90000;

    GridConfig gc;
    gc.num_cells = N_total;
    gc.dx        = dx;
    gc.courant   = Q;
    gc.c_speed   = c0;
    Grid grid(gc);

    std::vector<Layer> layers;
    double x_pos = struct_start * dx;

    for (int p = 0; p < n_periods; ++p) {
        layers.push_back({x_pos, x_pos + L_TiO2, eps_TiO2, 1.0, 0.0});
        x_pos += L_TiO2;
        layers.push_back({x_pos, x_pos + L_SiO2, eps_SiO2, 1.0, 0.0});
        x_pos += L_SiO2;
    }

    double cavity_x0 = x_pos;
    x_pos += L_cavity;
    double cavity_x1 = x_pos;

    for (int p = 0; p < n_periods; ++p) {
        layers.push_back({x_pos, x_pos + L_SiO2, eps_SiO2, 1.0, 0.0});
        x_pos += L_SiO2;
        layers.push_back({x_pos, x_pos + L_TiO2, eps_TiO2, 1.0, 0.0});
        x_pos += L_TiO2;
    }

    applyLayers(grid, layers);

    PMLConfig pml_cfg;
    pml_cfg.thickness = pml_N;
    pml_cfg.order     = 3;
    pml_cfg.delta     = 1e-8;
    applyPML(grid, pml_cfg);

    auto wf = makeGaussianPulseForRange(lam_min, lam_max, c0, 6.0);
    Source source{src_pos, wf, SourceMode::Soft};

    std::vector<int> snap_times = {3000, 8000, 15000, 30000, 50000, 80000};
    std::vector<std::vector<double>> snap_Ey;
    std::vector<int> snap_labels;

    std::vector<double> x_nm(N_total);
    for (int i = 0; i < N_total; ++i) x_nm[i] = grid.xE(i) * 1e9;

    int snap_idx = 0;
    for (int n = 0; n < num_steps; ++n) {
        FDTDSolver::updateE(grid);
        source.inject(grid, n);
        FDTDSolver::updateH(grid);

        if (snap_idx < (int)snap_times.size() && n == snap_times[snap_idx]) {
            std::vector<double> ey(N_total);
            for (int i = 0; i < N_total; ++i) ey[i] = grid.Ey[i];
            snap_Ey.push_back(ey);
            snap_labels.push_back(n);
            ++snap_idx;
        }
    }

    double max_abs = 1e-12;
    for (const auto& s : snap_Ey)
        for (double v : s)
            max_abs = std::max(max_abs, std::abs(v));

    using namespace matplot;
    auto fig = figure(true);
    fig->size(1200, 600);
    hold(on);

    for (size_t s = 0; s < snap_Ey.size(); ++s) {
        auto p = plot(x_nm, snap_Ey[s]);
        p->display_name("n=" + std::to_string(snap_labels[s]));
        p->line_width(1.2);
    }

    auto l1 = plot({cavity_x0 * 1e9, cavity_x0 * 1e9}, {-max_abs, max_abs});
    l1->line_style("--"); l1->color("gray"); l1->line_width(1.0);
    l1->display_name("Cavity edges");

    auto l2 = plot({cavity_x1 * 1e9, cavity_x1 * 1e9}, {-max_abs, max_abs});
    l2->line_style("--"); l2->color("gray"); l2->line_width(1.0);

    xlabel("x (nm)");
    ylabel("E_y");
    matplot::title("Gaussian pulse: field evolution in Bragg microcavity");
    matplot::legend();
    show();
}

int main() {
    std::cout << "=== Task 5: Bragg Microcavity ===\n\n";

    double lam_center = 650.0;          // nm
    double cavity_nm  = lam_center / 2; // half-wave air cavity

    {
        std::cout << "[1] R, T vs mirror periods\n";
        std::vector<int> n_per_list = {8};

        using namespace matplot;

        for (int np : n_per_list) {
            auto r = runBraggCavity(np, cavity_nm, lam_center,
                                    std::to_string(np) + " per/side");

            auto fig = figure(true);
            fig->size(1000, 600);
            hold(on);

            auto p1 = plot(r.lam_nm, r.R);
            p1->display_name("R"); p1->line_width(2);

            auto p2 = plot(r.lam_nm, r.T);
            p2->display_name("T"); p2->line_width(2);

            auto p3 = plot(r.lam_nm, r.RTsum);
            p3->display_name("R+T");
            p3->line_style("--");
            p3->line_width(1.2);
            p3->color("gray");

            xlabel("Wavelength (nm)");
            ylabel("Coefficient");
            matplot::title("Bragg cavity, " + r.label);
            matplot::legend();
            xlim({300.0, 900.0});
            ylim({0.0, 1.1});
            show();
        }
    }

    {
        std::cout << "\n[2] Cavity vs pure PhC\n";
        int n_per = 8;

        auto r_cav  = runBraggCavity(n_per, cavity_nm, lam_center, "Cavity 8+8");
        auto r_pure = runPurePhC(2 * n_per, lam_center);

        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);

        auto p1 = plot(r_cav.lam_nm, r_cav.T);
        p1->display_name("T: Bragg cavity");
        p1->line_width(2);

        auto p2 = plot(r_pure.lam_nm, r_pure.T);
        p2->display_name("T: Pure PhC");
        p2->line_width(2);
        p2->line_style("--");

        xlabel("Wavelength (nm)");
        ylabel("T");
        matplot::title("Transmission: Bragg cavity vs pure photonic crystal");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.1});
        show();
    }

    {
        std::cout << "\n[3] Gaussian field snapshots\n";
        runFieldSnapshotsGaussian(8, cavity_nm, lam_center);

        std::cout << "\n[4] CW field snapshots\n";
        runFieldSnapshotsCW(8, cavity_nm, lam_center);
    }

    std::cout << "\nTask 5 complete.\n";
}