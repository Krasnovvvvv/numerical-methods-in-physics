#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>

#include "Helpers/Plotter.h"
#include "Labs/Lab7/Tasks/TaskWaveBase.h"
#include "Labs/Lab7/Tasks/TaskWaveWeighted.h"
#include "Labs/Lab7/ExactSolutions/ExactWaveSolutions.h"
#include "Labs/Lab7/Analysis/WaveSpectrumAnalyzer.h"

TaskWaveBase::ExactFunc exact =
    [](double x, double t, const TaskWaveBase::PhysParams& phys) {
        return exact_piecewise_velocity(x, t, phys);
};

int main() {
    Plotter plotter;

    TaskWaveBase::PhysParams phys{
        .c     = 100.0,
        .L     = 1.0,
        .x0    = 0.3,
        .delta = 0.05,
        .v0    = 1.0
    };

    double tMax = 0.05;
    double h_x  = 0.01;

    std::cout << "\n========== TASK 3: Spectral analysis (FFT) ==========\n";

    double tau = 0.5 * h_x / phys.c;

    TaskWaveWeighted taskW(&plotter, exact,
                           TaskWaveBase::PLOT_SOLUTION, 0.25);
    taskW.run(phys, tMax, h_x, tau);

    WaveSpectrumAnalyzer analyzer(&plotter);

    int s_mid = static_cast<int>(taskW.getTime().size()) / 2;

    analyzer.analyzeAtTime(taskW, s_mid);
}
