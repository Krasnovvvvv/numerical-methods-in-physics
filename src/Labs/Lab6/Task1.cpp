#define _USE_MATH_DEFINES
#include "Labs/Lab1/SLAESolvers/ThomasSolver.h"
#include "Labs/Lab6/Tasks/TaskHeatEquation.h"
#include "Helpers/Plotter.h"
#include <cmath>

int main() {
    Plotter plotter;
    ThomasSolver solver;

    TaskHeatEquation task(solver, &plotter, 0.0);

    auto f   = [](double, double){ return 0.0; };
    auto u0  = [](double x){ return std::sin(M_PI * x); };
    auto nu1 = [](double){ return 0.0; };
    auto nu2 = [](double){ return 0.0; };
    auto exact = [](double x, double t){
        return std::exp(-M_PI*M_PI*t) * std::sin(M_PI*x);
    };

    double kappa = 1.0, l = 1.0, tMax = 1.0;
    int N = 1000, M = 1000;

    task.run(kappa, l, tMax, N, M, f, u0, nu1, nu2, exact);

}