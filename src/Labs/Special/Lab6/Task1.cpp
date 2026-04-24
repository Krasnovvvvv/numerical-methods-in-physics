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
    const int N = 800;
    const double dx = 1.0;
    const double Q  = 0.5;
    const int num_steps = 8000;

    const double f_min = 0.01;
    const double f_max = 0.5;
    const double f_center = 0.22;
    const double f_width  = 1.5;

    const int src_pos = pml_thickness + 30;
    const int mon_pos = N / 2;

    const int nf = 300;
    std::vector<double> freqs(nf);
    for (int i = 0; i < nf; ++i) {
        freqs[i] = f_min + (f_max - f_min) * i / (nf - 1.0);
    }

    auto wf = makeGaussianPulse(f_center, f_width, 6.0);

    auto E_total = runOnce(N, dx, Q, num_steps, src_pos, mon_pos, wf,
                           pml_thickness, pml_order, pml_delta, freqs);

    const int ref_pml_thick = 100;
    const int ref_pml_order = 3;
    const double ref_pml_delta = 1e-12;

    const int N_ref   = N + 2 * (ref_pml_thick - pml_thickness);
    const int src_ref = ref_pml_thick + 30;
    const int mon_ref = src_ref + (mon_pos - src_pos);

    auto E_inc = runOnce(N_ref, dx, Q, num_steps, src_ref, mon_ref, wf,
                         ref_pml_thick, ref_pml_order, ref_pml_delta, freqs);

    PMLResult result;
    result.freqs = freqs;
    result.R.assign(nf, 0.0);
    result.R_max = 0.0;

    double max_inc2 = 0.0;
    for (int k = 0; k < nf; ++k) {
        max_inc2 = std::max(max_inc2, std::norm(E_inc[k]));
    }

    const double threshold = max_inc2 * 1e-2;

    std::vector<double> R_raw(nf, 0.0);
    std::vector<bool> valid(nf, false);

    for (int k = 0; k < nf; ++k) {
        const double inc2 = std::norm(E_inc[k]);

        if (inc2 > threshold) {
            auto E_ref = E_total[k] - E_inc[k];
            R_raw[k] = std::min(std::norm(E_ref) / inc2, 1.0);
            valid[k] = true;
        }
    }

    const int half_window = 2;
    for (int k = 0; k < nf; ++k) {
        if (!valid[k]) {
            result.R[k] = 0.0;
            continue;
        }

        double sum = 0.0;
        int cnt = 0;
        const int k1 = std::max(0, k - half_window);
        const int k2 = std::min(nf - 1, k + half_window);

        for (int j = k1; j <= k2; ++j) {
            if (valid[j]) {
                sum += R_raw[j];
                ++cnt;
            }
        }

        result.R[k] = (cnt > 0) ? (sum / cnt) : 0.0;

        const bool in_band = (freqs[k] >= 0.03 && freqs[k] <= 0.30);
        if (in_band) {
            result.R_max = std::max(result.R_max, result.R[k]);
        }
    }

    return result;
}

int main() {
    std::cout << "=== Task 1: PML Reflection Study ===\n\n";

    const int fixed_width = 50;
    std::vector<int> orders = {0, 2, 3};
    std::vector<std::string> names = {
        "m=0 (const)",
        "m=2 (quadratic)",
        "m=3 (cubic)"
    };

    std::vector<PMLResult> results_rf;
    for (int ord : orders) {
        std::cout << "  R(f): order=" << ord
                  << ", width=" << fixed_width << "...\n";
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
        ylabel("R");
        matplot::title("PML Reflection vs Frequency (width=" +
                       std::to_string(fixed_width) + " cells)");
        matplot::legend();
        show();
    }

    std::vector<int> widths = {5, 10, 15, 20, 25, 30, 40, 50};
    std::vector<std::vector<double>> rmax_curves(orders.size());
    std::vector<double> widths_d(widths.begin(), widths.end());

    for (size_t j = 0; j < orders.size(); ++j) {
        for (int w : widths) {
            std::cout << "  R_max: order=" << orders[j]
                      << ", width=" << w << "...\n";
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
        ylabel("R_{max}");
        matplot::title("Maximum PML Reflection vs Width (working band, smoothed)");
        matplot::legend();
        show();
    }

    std::cout << "\nTask 1 complete.\n";
}