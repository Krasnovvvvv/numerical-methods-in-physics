#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
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

    TaskWaveBase::PhysParams phys{
        .c = 100.0,
        .L = 1.0,
        .x0 = 0.3,
        .delta = 0.05,
        .v0 = 1.0
    };

    double tMax = 0.05;
    double h_x = 0.01;

    std::cout << "\n========== TASK 1: Explicit Scheme (Warm-up) ==========\n";

    // Стабильный случай
    double tau_stable = 0.5 * h_x / phys.c;
    std::cout << "\nCase 1: Stable (lambda = 0.5, CFL satisfied)\n";

    TaskWaveExplicit taskStable(&plotter, exact, TaskWaveBase::PLOT_ERROR);
    taskStable.run(phys, tMax, h_x, tau_stable);

    // Нестабильный случай
    double tau_unstable = 2.0 * h_x / phys.c;
    std::cout << "\nCase 2: Unstable (lambda = 2.0, CFL violated)\n";

    TaskWaveExplicit taskUnstable(&plotter, nullptr);
    taskUnstable.run(phys, tMax, h_x, tau_unstable);

    return 0;
}
