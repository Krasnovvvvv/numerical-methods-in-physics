#include "Labs/Lab3/IntegralSolvers/AdaptiveSimpsonSolver.h"
#include "Labs/Lab3/Tasks/IntegrateTask.h"
#include <cmath>

int main() {
    constexpr double a = 0.0, b = 1.0, tol = 1e-13;

    // replace the variable by the formula: x = t / (1 - t); t ∈ [0, 1)
    auto integrand = [](double t) -> double {
        if (t <= 0.0) return 0.0;
        if (t >= 1.0) return 0.0;
        double x = t / (1 - t);
        double dxdt = 1.0 / ((1 - t) * (1 - t));
        return std::exp(-std::pow(x, 3)) * std::sin(x) * std::log(x) * dxdt;
    };

    AdaptiveSimpsonSolver adSimpSolver;

    Plotter plotter;

    IntegrateTask task_adSimpson(adSimpSolver, &plotter, 1);

    task_adSimpson.run(integrand, a, b, tol);


}