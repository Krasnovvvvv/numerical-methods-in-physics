/// Задача 0: Реализация 1D FDTD алгоритма
///
/// Требуемые графики:
/// 1. E_y(x) в разные моменты: Гаусс, soft source, БЕЗ PML — видны отражения от границ
/// 2. E_y(x) в разные моменты: Гаусс, soft source, С PML — волна поглощается
/// 3. CW-источник, soft source, с PML — стационарная волна
/// 4. Сравнение soft source vs current source — показать что разницы нет

#include "Labs/Special/Lab6/Base/Grid.h"
#include "Labs/Special/Lab6/Base/Source.h"
#include "Labs/Special/Lab6/Solver/FDTDSolver.h"
#include "Labs/Special/Lab6/Base/PML.h"

#include <matplot/matplot.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace fdtd;

struct Snapshot {
    int step;
    std::vector<double> Ey;
};

std::vector<Snapshot> runAndCollect(
    int N, double dx, double Q, int num_steps,
    int src_pos, Waveform wf, SourceMode mode,
    bool usePML, int pml_N, int snap_interval)
{
    GridConfig gc;
    gc.num_cells = N;
    gc.dx = dx;
    gc.courant = Q;
    // c_speed=1.0 по умолчанию (нормализованные единицы)
    Grid grid(gc);

    if (usePML) {
        PMLConfig pml;
        pml.thickness = pml_N;
        pml.order = 3;
        pml.delta = 1e-6;
        applyPML(grid, pml);
    }

    Source source{src_pos, wf, mode};
    std::vector<Snapshot> snapshots;

    for (int n = 0; n < num_steps; ++n) {
        FDTDSolver::updateE(grid);
        source.inject(grid, n);
        FDTDSolver::updateH(grid);

        if (n % snap_interval == 0 && n > 0) {
            Snapshot s;
            s.step = n;
            s.Ey.resize(N);
            for (int i = 0; i < N; ++i) s.Ey[i] = grid.Ey[i];
            snapshots.push_back(s);
        }
    }
    return snapshots;
}

void plotSnapshots(const std::vector<double>& x,
                   const std::vector<Snapshot>& snaps,
                   const std::string& title_str) {
    using namespace matplot;

    auto fig = figure(true);
    fig->size(1000, 500);
    hold(on);

    for (const auto& s : snaps) {
        auto p = plot(x, s.Ey);
        p->display_name("n = " + std::to_string(s.step));
        p->line_width(1.5);
    }
    xlabel("x (cells)");
    ylabel("E_y");
    matplot::title(title_str);
    matplot::legend()->location(legend::general_alignment::topleft);

    show();
}

int main() {
    std::cout << "=== Task 0: 1D FDTD in Free Space ===\n\n";

    int N = 500;
    double dx = 1.0;
    double Q = 0.5;
    int pml_N = 40;
    double freq = 0.05;

    std::vector<double> x(N);
    for (int i = 0; i < N; ++i) x[i] = i * dx;

    // ================================================================
    // 1. Гауссов импульс, soft source, БЕЗ PML — видны отражения
    // ================================================================
    {
        std::cout << "[1] Gauss pulse, soft source, NO PML\n";
        auto wf = makeGaussianPulse(freq, 0.03, 5.0);
        auto snaps = runAndCollect(N, dx, Q, 1200, 100, wf,
                                   SourceMode::Soft, false, 0, 300);
        plotSnapshots(x, snaps,
            "Gauss pulse, soft source, NO PML (reflections visible)");
    }

    // ================================================================
    // 2. Гауссов импульс, soft source, С PML — поглощение на границах
    // ================================================================
    {
        std::cout << "[2] Gauss pulse, soft source, WITH PML\n";
        auto wf = makeGaussianPulse(freq, 0.03, 5.0);
        auto snaps = runAndCollect(N, dx, Q, 1200, pml_N + 30, wf,
                                   SourceMode::Soft, true, pml_N, 300);
        plotSnapshots(x, snaps,
            "Gauss pulse, soft source, WITH PML (absorbed at boundaries)");
    }

    // ================================================================
    // 3. CW-источник, soft source, С PML — стационарная волна
    // ================================================================
    {
        std::cout << "[3] CW source, soft source, WITH PML\n";
        auto wf = makeCW(freq, 40.0, 3.0);
        auto snaps = runAndCollect(N, dx, Q, 800, pml_N + 30, wf,
                                   SourceMode::Soft, true, pml_N, 100);
        plotSnapshots(x, snaps,
            "CW source, soft source, WITH PML");
    }

    // ================================================================
    // 4. Сравнение soft vs current source — идентичны
    // ================================================================
    {
        std::cout << "[4] Comparison: soft vs current source\n";
        int src_pos = pml_N + 30;
        int snap_time = 400;
        auto wf = makeGaussianPulse(freq, 0.03, 5.0);

        auto snaps_soft = runAndCollect(N, dx, Q, 500, src_pos, wf,
                                        SourceMode::Soft, true, pml_N, snap_time);
        auto snaps_curr = runAndCollect(N, dx, Q, 500, src_pos, wf,
                                        SourceMode::Current, true, pml_N, snap_time);

        if (!snaps_soft.empty() && !snaps_curr.empty()) {
            using namespace matplot;

            auto fig = figure(true);
            fig->size(1000, 500);
            hold(on);

            auto p1 = plot(x, snaps_soft[0].Ey);
            p1->display_name("Soft source");
            p1->line_width(2.0);

            auto p2 = plot(x, snaps_curr[0].Ey);
            p2->display_name("Current source");
            p2->line_width(1.5);
            p2->line_style("--");

            xlabel("x (cells)");
            ylabel("E_y");
            matplot::title("Soft vs Current source at step " +
                           std::to_string(snaps_soft[0].step));
            matplot::legend();

            show();
        }
    }

    std::cout << "\nTask 0 complete.\n";
    std::cout << "\nВыводы:\n";
    std::cout << "  1. Без PML — волна отражается от границ расчётной области\n";
    std::cout << "  2. С PML — волна поглощается, отражений нет\n";
    std::cout << "  3. CW-источник формирует стационарную синусоидальную волну\n";
    std::cout << "  4. Soft source и current source дают идентичные результаты —\n";
    std::cout << "     оба прозрачны для отражённых волн\n";
}