#include <cmath>
#include <vector>
#include <memory>
#include <iostream>

// PDE solvers
#include "Labs/Special/Lab3/ConvectionSolvers/FTCSSolver.h"
#include "Labs/Special/Lab3/ConvectionSolvers/LaxWendroffSolver.h"
#include "Labs/Special/Lab3/ConvectionSolvers/RichtmyerSolver.h"
#include "Labs/Special/Lab3/ConvectionSolvers/MacCormackSolver.h"
#include "Labs/Special/Lab3/ConvectionSolvers/Upwind1Solver.h"
#include "Labs/Special/Lab3/ConvectionSolvers/Upwind2Solver.h"
#include "Labs/Special/Lab3/ConvectionSolvers/BTCSSolver.h"

// Task
#include "Labs/Special/Lab3/Tasks/ConvectionTask.h"

// Linear solver for BTCS
#include "Labs/Lab1/SLAESolvers/ThomasSolver.h"

int main() {
    using namespace special;

    const double xL = 0.0;
    const double xR = 1.0;
    const std::size_t Nx = 201;      // число узлов по x

    const double t0 = 0.0;
    const double tN = 10.0;

    std::vector<double> outputTimes{0.1, 0.5, 1.0, 5.0, 10.0};

    // скорость переноса (можно менять, но c = u*dt/dx будет подбираться ниже)
    const double u = 0.1;

    auto phi = [](double x) {
        return std::cos(M_PI * x / 2.0);
    };

    auto exact_solution = [=](double x, double t) {
        double arg = x - u * t;
        return std::cos(M_PI * arg / 2.0);
    };

    auto leftBC  = [=](double t) { return exact_solution(xL, t); };
    auto rightBC = [=](double t) { return exact_solution(xR, t); };

    Plotter plotter;

    ThomasSolver linearSolver;
    BTCSSolver btcsSolver(linearSolver);

    FTCSSolver        ftcsSolver;
    LaxWendroffSolver lwSolver;
    RichtmyerSolver   richtSolver;
    MacCormackSolver  macSolver;
    Upwind1Solver     up1Solver;
    Upwind2Solver     up2Solver;

    std::vector<IConvectionSolver*> solvers = {
        &ftcsSolver,
        &lwSolver,
        &richtSolver,
        &macSolver,
        &up1Solver,
        &up2Solver,
        &btcsSolver
    };

    std::vector<double> c_values{0.1, 0.5, 1.0};

    for (double c : c_values) {
        double dx = (xR - xL) / (Nx - 1);
        double dt = c * dx / u;

        std::cout << "\n=============================\n";
        std::cout << "Convection number c = " << c
                  << ", dt = " << dt << "\n";

        ConvectionTaskParams params;
        params.xL = xL;
        params.xR = xR;
        params.Nx = Nx;
        params.t0 = t0;
        params.tN = tN;
        params.dt = dt;
        params.convectionVelocity = u;
        params.u0    = phi;
        params.exact = exact_solution;
        params.leftBC  = leftBC;
        params.rightBC = rightBC;

        for (IConvectionSolver* s : solvers) {
            ConvectionTask task(*s, &plotter, "x", "U");
            task.run(params, outputTimes);
        }
    }
}