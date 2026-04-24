#include "Labs/Special/Lab6/Base/Grid.h"
#include "Labs/Special/Lab6/Base/Source.h"
#include "Labs/Special/Lab6/Base/PML.h"
#include "Labs/Special/Lab6/Base/Monitor.h"
#include "Labs/Special/Lab6/Base/Material.h"
#include "Labs/Special/Lab6/Simulation/Simulation.h"

#include <matplot/matplot.h>

#include <iostream>
#include <vector>
#include <cmath>

using namespace fdtd;

int main() {
    std::cout << "=== Task 2: Reflection from Infinite Dielectric ===\n\n";

    const double c0         = 3e8;
    const double lambda_min = 300e-9;
    const double lambda_max = 900e-9;
    const double n_quartz   = 1.45;
    const double eps_quartz = n_quartz * n_quartz;

    const double dx    = 10e-9;
    const int pml_N    = 40;
    const int N_domain = 400;
    const int N_total  = N_domain + 2 * pml_N;
    const double Q     = 0.5;
    const int num_steps = 30000;

    int src_pos       = pml_N + 50;
    int mon_ref       = pml_N + 100;
    int interface_pos = pml_N + N_domain / 2;
    int mon_trans     = interface_pos + 50;

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

    double x_interface = interface_pos * dx;
    double x_end       = N_total * dx;
    sim.addLayer({x_interface, x_end, eps_quartz, 1.0, 0.0});

    sim.setSource(src_pos, makeGaussianPulseForRange(lambda_min, lambda_max, c0, 6.0),
                  SourceMode::Soft);
    sim.addMonitor(mon_ref, freqs);
    sim.addMonitor(mon_trans, freqs);
    sim.build();

    std::cout << "  Running main simulation (" << num_steps << " steps)...\n";
    sim.run();

    std::cout << "  Running normalization...\n";
    auto norm_monitors = sim.runNormalization();
    auto rt = sim.computeRT(norm_monitors, c0);

    double R_theory = std::pow((1.0 - n_quartz) / (1.0 + n_quartz), 2);
    double T_theory = 1.0 - R_theory;

    std::vector<double> T_power(rt.T.size());
    for (size_t k = 0; k < rt.T.size(); ++k)
        T_power[k] = n_quartz * rt.T[k];

    using namespace matplot;
    std::vector<double> lam_nm(rt.wavelengths.size());
    for (size_t i = 0; i < lam_nm.size(); ++i)
        lam_nm[i] = rt.wavelengths[i] * 1e9;

    std::vector<double> R_th(lam_nm.size(), R_theory);
    std::vector<double> T_th(lam_nm.size(), T_theory);

    auto fig = figure(true);
    fig->size(1000, 600);
    hold(on);

    auto p1 = plot(lam_nm, rt.R);
    p1->display_name("FDTD R");
    p1->line_width(2);

    auto p2 = plot(lam_nm, R_th);
    p2->display_name("Fresnel R = " + std::to_string(R_theory).substr(0, 6));
    p2->line_style("--");
    p2->line_width(2);
    p2->color("red");

    auto p3 = plot(lam_nm, T_power);
    p3->display_name("FDTD T");
    p3->line_width(2);

    auto p4 = plot(lam_nm, T_th);
    p4->display_name("Fresnel T = " + std::to_string(T_theory).substr(0, 6));
    p4->line_style("--");
    p4->line_width(2);
    p4->color("green");

    xlabel("Wavelength (nm)");
    ylabel("Coefficient");
    matplot::title("Reflection & Transmission: SiO2 half-space (n=1.45)");
    matplot::legend();
    xlim({300.0, 900.0});
    ylim({0.0, 1.1});
    show();

    std::cout << "\nTask 2 complete.\n";
}