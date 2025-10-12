#include "Helpers/Plotter.h"
#include "Labs/Lab4/ODESolvers/EulerSolver.h"
#include "Labs/Lab4/ODESolvers/ImprovedEulerSolver.h"
#include "Labs/Lab4/ODESolvers/RungeKutta4Solver.h"
#include "Labs/Lab4/Tasks/ODETask.h"
#include "Labs/Lab4/Tasks/Strategies/ClassicStrategy.h"
#include "Labs/Lab4/Tasks/Strategies/SolverDependenceStrategy.h"
#include "Labs/Lab4/Tasks/Strategies/StepDependenceStrategy.h"
#include "Helpers/ODEFuncs.h"

int main() {

    constexpr double t0 = 0.0, tn = 10.0, h = 0.001, tol = 1e-8;
    std::vector<double> y0 = {1.0, 0.0};

    Plotter plotter;
    EulerSolver euler;
    ImprovedEulerSolver improvedEuler;
    RungeKutta4Solver rk4;

    // --- CLASSIC ---
    std::cout << "\n--- CLASSIC MODE ---\n";
        {
        auto strategy = std::make_unique<ClassicStrategy>(&rk4);
        ODETask task(std::move(strategy));
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 1, &plotter); // x(t)
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 2, &plotter); // v(t)
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 3, &plotter); // error
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 4, &plotter); // phase traj
    }

    // --- DEPENDENCE_STEP ---
    std::cout << "\n--- DEPENDENCE_STEP MODE ---\n";
    {
        std::vector<double> steps = {0.01, 0.001};
        auto strategy = std::make_unique<StepDependenceStrategy>(&rk4);
        ODETask task(std::move(strategy));
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 1, &plotter, steps);
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 2, &plotter, steps);
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 3, &plotter, steps);
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 4, &plotter, steps);
    }

    // --- DEPENDENCE_SOLVER ---
    std::cout << "\n--- DEPENDENCE_SOLVER MODE ---\n";
    {
        std::vector<IODESolver*> solvers = {&euler, &improvedEuler, &rk4};
        auto strategy = std::make_unique<SolverDependenceStrategy>();
        ODETask task(std::move(strategy));
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 1, &plotter, {}, {}, solvers);
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 2, &plotter, {}, {}, solvers);
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 3, &plotter, {}, {}, solvers);
        task.run(harmonicOscillator, y0, t0, tn, h, tol, 4, &plotter, {}, {}, solvers);
    }
}
