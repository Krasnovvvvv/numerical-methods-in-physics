/// Задача 1: Изучение отражения от PML
///
/// Метод: два прогона с ОДНИМ монитором слева от правого PML.
///   Прогон 1: PML слева и справа → E_total(f) содержит падающую + отражённую от PML
///   Прогон 2: PML слева, а справа — очень длинная область (или толстый идеальный PML)
///             → E_inc(f) чисто падающая волна
///
/// Но проще: используем Simulation с фиктивной структурой (ε=1).
/// Нормировка — прогон с "идеальным" PML (очень толстый, m=3).
/// Разница между E_total и E_inc даёт E_reflected от PML.
///
/// Ещё проще: один прогон. Монитор ставим ПОЗАДИ источника (между источником и левым PML).
/// Падающая волна идёт ВПРАВО → монитор её не видит в прямом направлении.
/// Отражённая от правого PML идёт ВЛЕВО → монитор её ловит.
/// Но soft source прозрачен, значит отражение проходит через него.
/// Однако E_inc тоже ловится при t=0, когда импульс рождается.
///
/// Самый чистый метод: два прогона, монитор между источником и правым PML.

#include "Labs/Special/Lab6/Base/Grid.h"
#include "Labs/Special/Lab6/Base/Source.h"
#include "Labs/Special/Lab6/Solver/FDTDSolver.h"
#include "Labs/Special/Lab6/Base/PML.h"
#include "Labs/Special/Lab6/Base/Monitor.h"

#include <matplot/matplot.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <string>

using namespace fdtd;

struct PMLResult {
    std::vector<double> freqs;
    std::vector<double> R;
    double R_max;
};

/// Один прогон FDTD с PML и DFT-монитором. Возвращает комплексный спектр E(f).
std::vector<std::complex<double>> runOnce(
    int N, double dx, double Q, int num_steps,
    int src_pos, int mon_pos, const Waveform& wf,
    int pml_thickness, int pml_order, double pml_delta,
    const std::vector<double>& freqs)
{
    GridConfig gc;
    gc.num_cells = N;
    gc.dx = dx;
    gc.courant = Q;
    Grid grid(gc);

    PMLConfig pml;
    pml.thickness = pml_thickness;
    pml.order = pml_order;
    pml.delta = pml_delta;
    applyPML(grid, pml);

    Source source{src_pos, wf, SourceMode::Soft};
    DFTMonitor monitor(mon_pos, freqs);

    for (int n = 0; n < num_steps; ++n) {
        FDTDSolver::updateE(grid);
        source.inject(grid, n);
        monitor.accumulateAfterE(grid, n);
        FDTDSolver::updateH(grid);
    }

    return monitor.E_dft();
}

PMLResult runPMLExperiment(int pml_thickness, int pml_order, double pml_delta = 1e-6) {
    int N = 800;
    double dx = 1.0;
    double Q  = 0.5;
    int num_steps = 4000;  // импульс короткий, 4000 шагов достаточно

    // Широкополосный импульс: покрывает 0.01 — 0.33 (основная часть
    // разрешимого диапазона, f_Nyquist = 0.5).
    double f_center = 0.22;
    double f_width  = 1.5;   // 1/τ — короткий импульс, широкий спектр

    int src_pos = pml_thickness + 30;
    int mon_pos = (N / 2);

    int nf = 300;
    std::vector<double> freqs(nf);
    double f_plot_max = 0.35;
    for (int i = 0; i < nf; ++i) {
        freqs[i] = 0.005 + f_plot_max * i / (nf - 1.0);
    }

    auto wf = makeGaussianPulse(f_center, f_width, 6.0);

    // Прогон 1: с тестируемым PML → E_total = E_inc + E_ref
    auto E_total = runOnce(N, dx, Q, num_steps, src_pos, mon_pos, wf,
                           pml_thickness, pml_order, pml_delta, freqs);

    // Прогон 2 (нормировка): с очень толстым идеальным PML → E_inc (отражение ~0)
    int ref_pml_thick = 100;  // толстый PML для эталона
    int ref_pml_order = 3;
    double ref_pml_delta = 1e-12;

    // Для нормировки нужна бóльшая область, чтобы поместился толстый PML
    int N_ref = N + 2 * (ref_pml_thick - pml_thickness);
    int src_ref = ref_pml_thick + 30;
    int mon_ref = src_ref + (mon_pos - src_pos);  // тот же сдвиг от источника

    auto E_inc = runOnce(N_ref, dx, Q, num_steps, src_ref, mon_ref, wf,
                         ref_pml_thick, ref_pml_order, ref_pml_delta, freqs);

    // R(f) = |E_total - E_inc|^2 / |E_inc|^2
    PMLResult result;
    result.freqs = freqs;
    result.R.resize(nf);
    result.R_max = 0;

    // Относительный порог: считаем R только там, где спектр значим
    double max_inc2 = 0;
    for (int k = 0; k < nf; ++k) {
        double v = std::norm(E_inc[k]);
        if (v > max_inc2) max_inc2 = v;
    }
    double threshold = max_inc2 * 1e-4;

    for (int k = 0; k < nf; ++k) {
        double inc2 = std::norm(E_inc[k]);
        if (inc2 > threshold) {
            auto E_ref = E_total[k] - E_inc[k];
            result.R[k] = std::min(std::norm(E_ref) / inc2, 1.0);
        } else {
            result.R[k] = 0;
        }
        if (result.R[k] > result.R_max) result.R_max = result.R[k];
    }

    return result;
}

int main() {
    std::cout << "=== Task 1: PML Reflection Study ===\n\n";

    // ================================================================
    // График 1: R(f) при фиксированной ширине для разных профилей
    // ================================================================
    int fixed_width = 20;
    std::vector<int> orders = {0, 2, 3};
    std::vector<std::string> names = {"m=0 (const)", "m=2 (quadratic)", "m=3 (cubic)"};

    std::vector<PMLResult> results_rf;
    for (int ord : orders) {
        std::cout << "  R(f): order=" << ord << ", width=" << fixed_width << "...\n";
        results_rf.push_back(runPMLExperiment(fixed_width, ord));
        std::cout << "    R_max = " << results_rf.back().R_max << "\n";
    }

    {
        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);
        for (size_t i = 0; i < results_rf.size(); ++i) {
            std::vector<double> R_plot(results_rf[i].R.size());
            for (size_t k = 0; k < R_plot.size(); ++k)
                R_plot[k] = std::max(results_rf[i].R[k], 1e-16);
            auto p = semilogy(results_rf[i].freqs, R_plot);
            p->display_name(names[i]);
            p->line_width(2);
        }
        xlabel("Frequency (normalized)");
        ylabel("|R|^2");
        matplot::title("PML Reflection vs Frequency (width=" +
                       std::to_string(fixed_width) + " cells)");
        matplot::legend();
        show();
    }

    // ================================================================
    // График 2: R_max от ширины PML для разных профилей
    // ================================================================
    std::vector<int> widths = {5, 10, 15, 20, 25, 30, 40, 50};
    std::vector<std::vector<double>> rmax_curves(orders.size());
    std::vector<double> widths_d(widths.begin(), widths.end());

    for (size_t j = 0; j < orders.size(); ++j) {
        for (int w : widths) {
            std::cout << "  R_max: order=" << orders[j] << ", width=" << w << "...\n";
            auto result = runPMLExperiment(w, orders[j]);
            rmax_curves[j].push_back(std::max(result.R_max, 1e-16));
        }
    }

    {
        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 600);
        hold(on);
        for (size_t i = 0; i < orders.size(); ++i) {
            auto p = semilogy(widths_d, rmax_curves[i]);
            p->display_name(names[i]);
            p->line_width(2);
            p->marker(line_spec::marker_style::circle);
        }
        xlabel("PML width (cells)");
        ylabel("Max |R|^2");
        matplot::title("Maximum PML Reflection vs Width");
        matplot::legend();
        show();
    }

    std::cout << "\nTask 1 complete.\n";
    std::cout << "\nВыводы:\n";
    std::cout << "  - Постоянный профиль (m=0) даёт наибольшее отражение\n";
    std::cout << "  - Кубический профиль (m=3) — наименьшее\n";
    std::cout << "  - Увеличение толщины PML снижает отражение\n";
}