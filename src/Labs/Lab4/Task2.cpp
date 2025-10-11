#include "Base/IODESolver.h"
#include "Labs/Lab4/ODESolvers/EulerSolver.h"
#include "Labs/Lab4/ODESolvers/ImprovedEulerSolver.h"
#include "Labs/Lab4/ODESolvers/RungeKutta4Solver.h"
#include "Helpers/Plotter.h"
#include "Labs/Lab4/TAsks/ODETask.h"
#include <vector>

// task parameters
constexpr double b = 1.0;
constexpr double omega = 10.0;

std::vector<double> euler_rhs(double t, const std::vector<double>& y) {
    return {
        y[1],
        - (2 * b / t) * y[1] - omega * omega * y[0]
    };
}

int main() {

    double t0 = 1.0;
    double tn = 100.0;
    double tol = 1e-8;

    double tau = 1e-3;

    // x(1) = 1, v(1) = 1
    std::vector<double> y0 = {1.0, 1.0};
    Plotter plotter;

    EulerSolver euler;
    ImprovedEulerSolver improvedEuler;
    RungeKutta4Solver rk4;
    ODETask taskEuler(euler, &plotter);
    ODETask taskImprovedEuler(improvedEuler, &plotter);
    ODETask taskRK4(rk4, &plotter);

    taskRK4.run(euler_rhs, y0, t0, tn, tau);

}
