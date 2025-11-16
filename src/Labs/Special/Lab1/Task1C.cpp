#include "Labs/Special/Lab1/ODESolvers/GearBDFSolver.h"
#include "Labs/Lab4/Tasks/Strategies/StepDependenceStrategy.h"
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
    GearBDFSolver solver(1);

    // --- Solution by implicit Euler (Geer-1) solver ---
    std::cout << "\n--- Solution by implicit Euler solver ---\n";
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
        tn = 0.1;
        StepDependenceStrategy strategy(&solver, {"y", "z"});
        std::vector<double> steps = {0.01, 0.02, 0.025};

        // --- solutions for different h ---
        strategy.run(rhs, y0, t0, tn, h, tol, 1, &plotter, steps, {}, {});

        // --- |y_exact - y_num| for different h
        strategy.run(rhs, y0, t0, tn, h, tol, 3, &plotter, steps, {}, {});
    }
}