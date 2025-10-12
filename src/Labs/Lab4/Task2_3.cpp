#include "Helpers/Plotter.h"
#include "Labs/Lab4/ODESolvers/EulerSolver.h"
#include "Labs/Lab4/ODESolvers/ImprovedEulerSolver.h"
#include "Labs/Lab4/ODESolvers/RungeKutta4Solver.h"
#include "Labs/Lab4/Tasks/ODETask.h"
#include "Labs/Lab4/Tasks/Strategies/ClassicStrategy.h"
#include "Labs/Lab4/Tasks/Strategies/SolverDependenceStrategy.h"
#include "Labs/Lab4/Tasks/Strategies/StepDependenceStrategy.h"
#include "Labs/Lab4/Tasks/Strategies/ToleranceDependenceStrategy.h"
#include "Labs/Lab4/Tasks/Strategies/CompareWithExpectedStrategy.h"
#include "Helpers/ODEFuncs.h"

inline std::pair<double, double> getC1C2(double t0, double x0, double v0, double b, double omega) {
    double tb = std::pow(t0, -b);
    double coswt = std::cos(omega * t0);
    double sinwt = std::sin(omega * t0);

    double rhs1 = x0 / tb;
    double rhs2 = v0 / tb;

    double M11 = coswt;
    double M12 = sinwt;
    double M21 = -b / t0 * coswt - omega * sinwt;
    double M22 = -b / t0 * sinwt + omega * coswt;
    double det = M11 * M22 - M12 * M21;
    double C1 = (rhs1 * M22 - rhs2 * M12) / det;
    double C2 = (rhs2 * M11 - rhs1 * M21) / det;
    return {C1, C2};
}

inline std::function<double(double)> make_analytic_solution(
    double t0, double x0, double v0, double b, double omega)
{
    auto coeffs = getC1C2(t0, x0, v0, b, omega);
    double C1 = coeffs.first, C2 = coeffs.second;
    return [=](double t) {
        return std::pow(t, -b) * (C1 * std::cos(omega * t) + C2 * std::sin(omega * t));
    };
}

int main() {

    constexpr double t0 = 1.0, tn = 100.0, h = 0.01, tol = 1e-8;
    std::vector<double> y0 = {1.0, 1.0};

    auto analytic = make_analytic_solution(t0, y0[0], y0[1], 1.0,10.0);

    Plotter plotter;
    EulerSolver euler;
    ImprovedEulerSolver improvedEuler;
    RungeKutta4Solver rk4;

    // --- Comparing with analytic solution ---
    {
        auto strategy = std::make_unique<CompareWithExpectedStrategy>(&rk4, analytic);
        ODETask task(std::move(strategy));

        task.run(oscillator, y0, t0, tn, h, tol, 1, &plotter); // comparing
        task.run(oscillator, y0, t0, tn, h, tol, 2, &plotter); // absolute error
    }

    // --- CLASSIC ---
    std::cout << "\n--- CLASSIC MODE ---\n";
        {
        auto strategy = std::make_unique<ClassicStrategy>(&rk4);
        ODETask task(std::move(strategy));
        task.run(oscillator, y0, t0, tn, h, tol, 1, &plotter); // x(t)
        task.run(oscillator, y0, t0, tn, h, tol, 2, &plotter); // v(t)
        task.run(oscillator, y0, t0, tn, h, tol, 3, &plotter); // error
        task.run(oscillator, y0, t0, tn, h, tol, 4, &plotter); // phase traj
    }

    // --- DEPENDENCE_STEP ---
    std::cout << "\n--- DEPENDENCE_STEP MODE ---\n";
    {
        std::vector<double> steps = {0.01, 0.001};
        auto strategy = std::make_unique<StepDependenceStrategy>(&rk4);
        ODETask task(std::move(strategy));
        task.run(oscillator, y0, t0, tn, h, tol, 1, &plotter, steps);
        task.run(oscillator, y0, t0, tn, h, tol, 2, &plotter, steps);
        task.run(oscillator, y0, t0, tn, h, tol, 3, &plotter, steps);
        task.run(oscillator, y0, t0, tn, h, tol, 4, &plotter, steps);
    }

    // --- DEPENDENCE_TOLERANCE ---
    std::cout << "\n--- DEPENDENCE_TOLERANCE MODE ---\n";
    {
        std::vector<double> tolerances = {1e-4, 1e-12};
        auto strategy = std::make_unique<ToleranceDependenceStrategy>(&euler);
        ODETask task(std::move(strategy));
        task.run(oscillator, y0, t0, tn, h, tol, 1, &plotter, {}, tolerances);
        task.run(oscillator, y0, t0, tn, h, tol, 2, &plotter, {}, tolerances);
        task.run(oscillator, y0, t0, tn, h, tol, 3, &plotter, {}, tolerances);
        task.run(oscillator, y0, t0, tn, h, tol, 4, &plotter, {}, tolerances);
    }

    // --- DEPENDENCE_SOLVER ---
    std::cout << "\n--- DEPENDENCE_SOLVER MODE ---\n";
    {
        std::vector<IODESolver*> solvers = {&euler, &improvedEuler, &rk4};
        auto strategy = std::make_unique<SolverDependenceStrategy>();
        ODETask task(std::move(strategy));
        task.run(oscillator, y0, t0, tn, h, tol, 1, &plotter, {}, {}, solvers);
        task.run(oscillator, y0, t0, tn, h, tol, 2, &plotter, {}, {}, solvers);
        task.run(oscillator, y0, t0, tn, h, tol, 3, &plotter, {}, {}, solvers);
        task.run(oscillator, y0, t0, tn, h, tol, 4, &plotter, {}, {}, solvers);
    }
}
