/// Задача 5: Брэгговская микрополость
///
/// Графики:
/// 1. R(λ), T(λ) микрополости — узкий пик пропускания внутри запрещённой зоны
/// 2. Сравнение с чисто периодической структурой (без дефекта)
/// 3. Снимки E_y(x) — распределение поля внутри микрополости

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

/// Результат
struct CavityResult {
    std::vector<double> lam_nm;
    std::vector<double> R, T;
    std::string label;
};

/// Прогон брэгговской микрополости:
///   [зеркало_L] | [полость (воздух)] | [зеркало_R]
/// Зеркало = N периодов (TiO₂/SiO₂), четвертьволновая укладка при λ_center
CavityResult runBraggCavity(int n_periods, double cavity_nm, double lambda_center_nm,
                             const std::string& label) {
    std::cout << "  Bragg: " << n_periods << " per/side, cavity="
              << cavity_nm << " nm, λ_c=" << lambda_center_nm << " nm\n";

    const double c0      = 3e8;
    const double lam_min = 300e-9;
    const double lam_max = 900e-9;
    const double n_SiO2  = 1.45;
    const double n_TiO2  = 2.28;
    const double eps_SiO2 = n_SiO2 * n_SiO2;
    const double eps_TiO2 = n_TiO2 * n_TiO2;

    double lam_c = lambda_center_nm * 1e-9;
    double L_SiO2 = lam_c / (4.0 * n_SiO2);
    double L_TiO2 = lam_c / (4.0 * n_TiO2);
    double L_cavity = cavity_nm * 1e-9;
    double period = L_SiO2 + L_TiO2;

    const double dx  = 5e-9;
    const int pml_N  = 40;
    const double Q   = 0.5;

    double struct_len = 2.0 * n_periods * period + L_cavity;
    int struct_cells  = (int)(struct_len / dx) + 10;
    int margin = 80;
    int N_total = 2 * pml_N + 2 * margin + struct_cells;

    int src_pos      = pml_N + 60;
    int struct_start = pml_N + margin;
    int struct_end   = struct_start + struct_cells;
    int mon_ref      = src_pos + 30;
    int mon_trans    = struct_end + 30;

    if (mon_trans + pml_N + 10 > N_total)
        N_total = mon_trans + pml_N + 20;

    int num_steps = 80000; // больше шагов для узкого резонанса
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

    // Левое зеркало: TiO₂ / SiO₂
    double x_pos = struct_start * dx;
    for (int p = 0; p < n_periods; ++p) {
        sim.addLayer({x_pos, x_pos + L_TiO2, eps_TiO2, 1.0, 0.0});
        x_pos += L_TiO2;
        sim.addLayer({x_pos, x_pos + L_SiO2, eps_SiO2, 1.0, 0.0});
        x_pos += L_SiO2;
    }

    // Полость (воздух) — пропускаем, ε=1 по умолчанию
    x_pos += L_cavity;

    // Правое зеркало: SiO₂ / TiO₂
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
    for (size_t k = 0; k < freqs.size(); ++k)
        res.lam_nm[k] = rt.wavelengths[k] * 1e9;

    return res;
}

/// Прогон чистого фотонного кристалла (без дефекта) для сравнения
CavityResult runPurePhC(int total_periods, double lambda_center_nm) {
    std::cout << "  Pure PhC: " << total_periods << " periods\n";

    const double c0      = 3e8;
    const double lam_min = 300e-9;
    const double lam_max = 900e-9;
    const double n_SiO2  = 1.45;
    const double n_TiO2  = 2.28;
    const double eps_SiO2 = n_SiO2 * n_SiO2;
    const double eps_TiO2 = n_TiO2 * n_TiO2;

    double lam_c = lambda_center_nm * 1e-9;
    double L_SiO2 = lam_c / (4.0 * n_SiO2);
    double L_TiO2 = lam_c / (4.0 * n_TiO2);
    double period = L_SiO2 + L_TiO2;

    const double dx  = 5e-9;
    const int pml_N  = 40;
    const double Q   = 0.5;

    double struct_len = total_periods * period;
    int struct_cells  = (int)(struct_len / dx) + 10;
    int margin = 80;
    int N_total = 2 * pml_N + 2 * margin + struct_cells;

    int src_pos      = pml_N + 60;
    int struct_start = pml_N + margin;
    int struct_end   = struct_start + struct_cells;
    int mon_ref      = src_pos + 30;
    int mon_trans    = struct_end + 30;

    if (mon_trans + pml_N + 10 > N_total)
        N_total = mon_trans + pml_N + 20;

    auto freqs = freqsForRange(lam_min, lam_max, 800, c0);

    Simulation::Config cfg;
    cfg.num_cells     = N_total;
    cfg.dx            = dx;
    cfg.courant       = Q;
    cfg.c_speed       = c0;
    cfg.num_steps     = 60000;
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
    for (size_t k = 0; k < freqs.size(); ++k)
        res.lam_nm[k] = rt.wavelengths[k] * 1e9;
    return res;
}

/// Прогон со снимками поля E_y(x)
void runFieldSnapshots(int n_periods, double cavity_nm, double lambda_center_nm) {
    std::cout << "  Field snapshots: " << n_periods << " per/side\n";

    const double c0      = 3e8;
    const double lam_min = 300e-9;
    const double lam_max = 900e-9;
    const double n_SiO2  = 1.45;
    const double n_TiO2  = 2.28;
    const double eps_SiO2 = n_SiO2 * n_SiO2;
    const double eps_TiO2 = n_TiO2 * n_TiO2;

    double lam_c   = lambda_center_nm * 1e-9;
    double L_SiO2  = lam_c / (4.0 * n_SiO2);
    double L_TiO2  = lam_c / (4.0 * n_TiO2);
    double L_cavity = cavity_nm * 1e-9;
    double period   = L_SiO2 + L_TiO2;

    const double dx  = 5e-9;
    const int pml_N  = 40;
    const double Q   = 0.5;

    double struct_len = 2.0 * n_periods * period + L_cavity;
    int struct_cells  = (int)(struct_len / dx) + 10;
    int margin = 80;
    int N_total = 2 * pml_N + 2 * margin + struct_cells;

    int src_pos      = pml_N + 60;
    int struct_start = pml_N + margin;
    int struct_end   = struct_start + struct_cells;

    if (struct_end + margin + pml_N > N_total)
        N_total = struct_end + margin + pml_N;

    int num_steps = 60000;

    GridConfig gc;
    gc.num_cells = N_total;
    gc.dx = dx;
    gc.courant = Q;
    gc.c_speed = c0;
    Grid grid(gc);

    // Левое зеркало
    double x_pos = struct_start * dx;
    std::vector<Layer> layers;
    for (int p = 0; p < n_periods; ++p) {
        layers.push_back({x_pos, x_pos + L_TiO2, eps_TiO2, 1.0, 0.0});
        x_pos += L_TiO2;
        layers.push_back({x_pos, x_pos + L_SiO2, eps_SiO2, 1.0, 0.0});
        x_pos += L_SiO2;
    }
    double cavity_center = x_pos + L_cavity / 2.0;
    x_pos += L_cavity;
    for (int p = 0; p < n_periods; ++p) {
        layers.push_back({x_pos, x_pos + L_SiO2, eps_SiO2, 1.0, 0.0});
        x_pos += L_SiO2;
        layers.push_back({x_pos, x_pos + L_TiO2, eps_TiO2, 1.0, 0.0});
        x_pos += L_TiO2;
    }
    applyLayers(grid, layers);

    PMLConfig pml_cfg;
    pml_cfg.thickness = pml_N;
    pml_cfg.order = 3;
    pml_cfg.delta = 1e-8;
    applyPML(grid, pml_cfg);

    auto wf = makeGaussianPulseForRange(lam_min, lam_max, c0, 6.0);
    Source source{src_pos, wf, SourceMode::Soft};

    // Собираем снимки в характерные моменты
    std::vector<int> snap_times = {2000, 5000, 10000, 20000, 40000, 55000};
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
            snap_idx++;
        }
    }

    // Показываем снимки поля
    using namespace matplot;
    auto fig = figure(true);
    fig->size(1200, 600);
    hold(on);

    for (size_t s = 0; s < snap_Ey.size(); ++s) {
        auto p = plot(x_nm, snap_Ey[s]);
        p->display_name("n=" + std::to_string(snap_labels[s]));
        p->line_width(1.2);
    }

    // Пометим положение полости
    double cav_nm = cavity_center * 1e9;
    auto vl = plot({cav_nm, cav_nm}, {-1.0, 1.0});
    vl->line_style("--"); vl->color("gray"); vl->line_width(1);
    vl->display_name("Cavity center");

    xlabel("x (nm)");
    ylabel("E_y");
    matplot::title("Field evolution in Bragg microcavity");
    matplot::legend();
    show();
}

int main() {
    std::cout << "=== Task 5: Bragg Microcavity ===\n\n";

    double lam_center = 550.0; // нм — центральная длина волны
    // Полость: полуволна воздуха = λ/2 = 275 нм
    double cavity_nm = lam_center / 2.0;

    // ================================================================
    // График 1: R(λ), T(λ) для разного числа зеркальных периодов
    // ================================================================
    {
        std::cout << "[1] R, T vs periods\n";
        std::vector<int> n_per_list = {3, 5, 8, 12};

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

            xlabel("Wavelength (nm)");
            ylabel("Coefficient");
            matplot::title("Bragg cavity, " + r.label + ", cavity=" +
                           std::to_string((int)cavity_nm) + " nm");
            matplot::legend();
            xlim({300.0, 900.0});
            ylim({0.0, 1.1});
            show();
        }
    }

    // ================================================================
    // График 2: Сравнение с чисто периодической структурой
    // ================================================================
    {
        std::cout << "\n[2] Cavity vs pure PhC\n";
        int n_per = 8;
        auto r_cav  = runBraggCavity(n_per, cavity_nm, lam_center,
                                      "Cavity 8+8");
        auto r_pure = runPurePhC(2 * n_per, lam_center); // 16 периодов

        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);

        auto p1 = plot(r_cav.lam_nm, r_cav.T);
        p1->display_name("T: Bragg cavity"); p1->line_width(2);
        auto p2 = plot(r_pure.lam_nm, r_pure.T);
        p2->display_name("T: Pure PhC (16 per)"); p2->line_width(2); p2->line_style("--");

        xlabel("Wavelength (nm)");
        ylabel("T");
        matplot::title("Transmission: Bragg cavity vs pure photonic crystal");
        matplot::legend();
        xlim({300.0, 900.0});
        ylim({0.0, 1.1});
        show();
    }

    // ================================================================
    // График 3: Снимки E_y(x) — эволюция поля
    // ================================================================
    {
        std::cout << "\n[3] Field snapshots\n";
        runFieldSnapshots(8, cavity_nm, lam_center);
    }

    std::cout << "\nTask 5 complete.\n";
    std::cout << "\nВыводы:\n";
    std::cout << "  - Внутри запрещённой зоны появляется узкий пик пропускания (мода полости)\n";
    std::cout << "  - Больше зеркальных периодов — уже резонанс (выше Q-фактор)\n";
    std::cout << "  - Чисто периодическая структура не имеет пика — только запрещённая зона\n";
    std::cout << "  - Поле локализуется внутри полости на резонансной частоте\n";
}