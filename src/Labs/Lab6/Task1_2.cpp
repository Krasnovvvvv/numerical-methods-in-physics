#define _USE_MATH_DEFINES
#include "Labs/Lab1/SLAESolvers/ThomasSolver.h"
#include "Labs/Lab6/Tasks/TaskHeatEquation.h"
#include "Helpers/Plotter.h"
#include <cmath>

int main() {
    Plotter plotter;
    ThomasSolver solver;

    int plotNumber = 2; // 1 — распределение, 2 — исследование ошибки
    double sigma = 1.0;

    TaskHeatEquation task(
        solver,
        &plotter,
        sigma,
        plotNumber == 1
            ? TaskHeatEquation::PlotMode::Distribution
            : TaskHeatEquation::PlotMode::ErrorStudy
    );

    auto f   = [](double, double){ return 0.0; };
    auto u0  = [](double x){ return std::sin(M_PI * x); };
    auto nu1 = [](double){ return 0.0; };
    auto nu2 = [](double){ return 0.0; };
    auto exact = [](double x, double t){
        return std::exp(-M_PI * M_PI * t) * std::sin(M_PI * x);
    };

    double kappa = 1.0, l = 1.0, tMax = 1.0;
    double h0 = 0.01;   // явный шаг по x
    double tau0 = 1.0 / 20; // явный шаг по t

    if (plotNumber == 1) {
        task.run(kappa, l, tMax, h0, tau0, f, u0, nu1, nu2, exact);
    } else {
        task.runErrorStudy(kappa, l, tMax, h0, tau0, f, u0, nu1, nu2, exact);
    }
}
