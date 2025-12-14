// Lab7_DenseWarmup.cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include "Helpers/Plotter.h"
#include "Labs/Lab7/Tasks/TaskWaveBase.h"
#include "Labs/Lab7/Tasks/TaskWaveExplicit.h"
#include "Labs/Lab7/ExactSolutions/ExactWaveSolutions.h"

TaskWaveBase::ExactFunc exact =
    [](double x, double t, const TaskWaveBase::PhysParams& phys) {
        return exact_piecewise_velocity(x, t, phys);
};

int main() {
    Plotter plotter;

    TaskWaveBase::PhysParams phys {
        .c     = 100.0,
        .L     = 1.0,
        .x0    = 0.3,
        .delta = 0.05,
        .v0    = 1.0
    };

    double tMax = 0.05;              // с
    double h_x  = 0.01;              // м
    double tau_stable   = 0.5 * h_x / phys.c;  // CFL
    double tau_unstable = 2.0 * h_x / phys.c;  // нарушение CFL

    // стабильный расчёт с точным решением
    TaskWaveExplicit taskStable(&plotter, exact);
    taskStable.run(phys, tMax, h_x, tau_stable);

    // прогонка с нарушением условия Куранта (без точного решения)
    TaskWaveExplicit taskUnstable(&plotter, nullptr);
    taskUnstable.run(phys, tMax, h_x, tau_unstable);

    std::cout << "Stable tau = "   << tau_stable
              << ", unstable tau = " << tau_unstable << "\n";

    return 0;
}
