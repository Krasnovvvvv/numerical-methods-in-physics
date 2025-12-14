// Lab7_LightHot.cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include "Helpers/Plotter.h"
#include "Base/ISolver.h"
#include "Labs/Lab1/SLAESolvers/ThomasSolver.h"

#include "Labs/Lab7/Tasks/TaskWaveBase.h"
#include "Labs/Lab7/Tasks/TaskWaveExplicit.h"
#include "Labs/Lab7/Tasks/TaskWaveTheta.h"
#include "Labs/Lab7/ExactSolutions/ExactWaveSolutions.h"

TaskWaveBase::ExactFunc exact =
    [](double x, double t, const TaskWaveBase::PhysParams& phys) {
        return exact_piecewise_velocity(x, t, phys);
};

int main() {
    Plotter     plotter;
    ThomasSolver solver;

    // физические параметры те же, что в "плотной разминке"
    TaskWaveBase::PhysParams phys {
        .c     = 100.0,
        .L     = 1.0,
        .x0    = 0.3,
        .delta = 0.05,
        .v0    = 1.0
    };

    double tMax = 0.05;   // с
    double h_x  = 0.01;   // м
    double tau  = 0.5 * h_x / phys.c;  // CFL для явной схемы

    // 1. Явная схема (для сравнения точности)
    TaskWaveExplicit explicitTask(&plotter, exact);
    explicitTask.run(phys, tMax, h_x, tau);

    // 2. Схема с весами: несколько значений sigma
    //    sigma < 1/4 — условная устойчивость, sigma >= 1/4 — безусловная
    std::vector<double> sigmas = {0.0, 0.1, 0.25, 0.4};

    for (double sigma : sigmas) {
        std::cout << "\n=== Theta-scheme run, sigma = " << sigma << " ===\n";
        TaskWaveTheta thetaTask(solver, &plotter, sigma, exact);
        thetaTask.run(phys, tMax, h_x, tau);
    }

    // 3. Отдельный эксперимент по устойчивости:
    //    увеличиваем шаг по времени при фиксированном h_x и фиксированной sigma.
    double sigma_test = 0.1;   // условно устойчивая зона (0 <= sigma < 1/4)
    std::cout << "\n=== Stability experiment for sigma = " << sigma_test << " ===\n";

    std::vector<double> taus = {
        0.25 * h_x / phys.c,   // сильно внутри CFL
        0.5  * h_x / phys.c,   // граница CFL явной схемы
        1.0  * h_x / phys.c,   // больше CFL
        2.0  * h_x / phys.c    // сильно больше
    };

    for (double tau_test : taus) {
        std::cout << "\n-- tau = " << tau_test << " --\n";
        TaskWaveTheta thetaStab(solver, &plotter,
                                sigma_test, exact);
        thetaStab.run(phys, tMax, h_x, tau_test);
    }
}
