#include "Base/IODESolver.h"
#include "Labs/Lab4/ODESolvers/EulerSolver.h"
#include "Labs/Lab4/ODESolvers/ImprovedEulerSolver.h"
#include "Labs/Lab4/ODESolvers/RungeKutta4Solver.h"
#include "Helpers/Plotter.h"
#include "Labs/Lab4/Tasks/OdeTask.h"
#include <vector>

// task parameters
constexpr double omega = 1.0;
std::vector<double> oscillator_rhs(double /*t*/, const std::vector<double>& y) {
    return {y[1], -omega * omega * y[0]};
}

int main() {
    double t0 = 0.0;
    double tn = 20.0;
    double h = 1e-3;
    std::vector<double> y0 = {1.0, 0.0}; // x(0)=1, v(0)=0

    Plotter plotter;

    EulerSolver euler;
    ImprovedEulerSolver improvedEuler;
    RungeKutta4Solver rk4;

    ODETask taskEuler(euler, &plotter);
    ODETask taskImprovedEuler(improvedEuler, &plotter);
    ODETask taskRK4(rk4, &plotter);

    taskRK4.run(oscillator_rhs, y0, t0, tn, h);

}
