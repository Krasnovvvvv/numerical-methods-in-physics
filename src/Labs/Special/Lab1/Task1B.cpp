#include "Labs/Lab4/ODESolvers/EulerSolver.h"
#include "Labs/Lab4/Tasks/Strategies/StepDependenceStrategy.h"
#include "Helpers/Plotter.h"
#include "Labs/Special/Lab1/Tasks/ODETask.h"


// --- Exact solutions ---
auto y_exact = [](double t) {
    return 4 * std::exp(-t) - 3 * std::exp(-1000 * t);
};

auto z_exact = [](double t) {
    return -2 * std::exp(-t) + 3 * std::exp(-1000 * t);
};

// --- System ---
auto rhs = [](double t, const std::vector<double>& Y) {
    double y = Y[0], z = Y[1];
    return std::vector<double>{
        998*y + 1998*z,
        -999*y - 1999*z
    };
};

int main() {
    std::vector<double> y0 = {1.0, 1.0}; // y(0) = 1, z(0) = 1
    double t0 = 0.0, tn = 0.1, h = 0.0001, tol = 1e-8;

    Plotter plotter;
    EulerSolver solver;

    // --- Solution by Euler solver ---
    std::cout << "\n--- Solution by Euler solver ---\n";
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
        tn = 0.01;
        StepDependenceStrategy strategy(&solver, {"y", "z"});
        std::vector<double> steps = {0.0005, 0.001, 0.002, 0.0025};

        // --- solutions for different h ---
        strategy.run(rhs, y0, t0, tn, h, tol, 1, &plotter, steps, {}, {});

        // --- |y_exact - y_num| for different h
        strategy.run(rhs, y0, t0, tn, h, tol, 3, &plotter, steps, {}, {});
    }
}