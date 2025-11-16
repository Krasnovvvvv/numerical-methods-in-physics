#include "Labs/Special/Lab1/ODESolvers/GearBDFSolver.h"
#include "Labs/Lab4/ODESolvers/RungeKutta4Solver.h"
#include "Labs/Lab4/Tasks/Strategies/SolverDependenceStrategy.h"
#include "Labs/Special/Lab1/Tasks/ODETask.h"
#include "Labs/Special/Lab1/ODEs/ODEs.h"
#include "Helpers/Plotter.h"

// --- Exact solutions ---
auto y_exact = rigidSystem::rigidSystemYSolution;

auto z_exact = rigidSystem::rigidSystemZSolution;

// --- System ---
auto rhs = rigidSystem::rigidSystem;

int main() {
    std::vector<double> y0 = {1.0, 1.0}; // y(0) = 1, z(0) = 1
    double t0 = 0.0, tn = 0.1, h = 0.001, tol = 1e-8;

    Plotter plotter;
    RungeKutta4Solver starter;
    GearBDFSolver solver(4, &starter);

    // --- Solution by Geer-4 solver with RK4 start ---
    std::cout << "\n--- Solution by Geer-4 solver with RK4 start ---\n";
    {
        // --- numeric vs analytic ---
        special::ODETask task1(solver, &plotter, 1, {"y", "z"});
        task1.run(rhs, y0, t0, tn, h, {y_exact, z_exact});

        // --- |y_exact - y_num| ---
        special::ODETask task2(solver, &plotter, 2, {"y", "z"});
        task2.run(rhs, y0, t0, tn, h, {y_exact, z_exact});
    }

    // --- Investigation of solution behavior ---
    {
        SolverDependenceStrategy strategy({"y", "z"});

        // --- Starters (Gear 1-3) ---
        GearBDFSolver starter1(1, &starter);
        GearBDFSolver starter2(2, &starter);
        GearBDFSolver starter3(3, &starter);

        // --- Solvers with different starters ---
        GearBDFSolver solver1(1, &starter1);
        GearBDFSolver solver2(2, &starter2);
        GearBDFSolver solver3(3, &starter3);

        std::vector<IODESolver*> solvers = {&solver1, &solver2, &solver3, &solver};

        // --- solutions for different solver ---
        strategy.run(rhs, y0, t0, tn, h, tol, 1, &plotter, {}, {}, solvers);

        // --- |y_exact - y_num| for different solver
        strategy.run(rhs, y0, t0, tn, h, tol, 3, &plotter, {}, {}, solvers);
    }
}