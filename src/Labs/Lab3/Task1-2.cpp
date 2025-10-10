#include "Labs/Lab3/IntegralSolvers/GaussSolver.h"
#include "Labs/Lab3/IntegralSolvers/CentralRectSolver.h"
#include "Labs/Lab3/IntegralSolvers/TrapezoidSolver.h"
#include "Labs/Lab3/IntegralSolvers/SimpsonSolver.h"
#include "Labs/Lab3/Tasks/IntegrateTask.h"
#include "Helpers/Plotter.h"
#include <cmath>

int main() {
    constexpr double a = 0.1, b = 0.4, tol = 1e-6;

    auto integrand = [](double x) {
        return std::atan(0.3 * std::pow(x, 4) - 4.0 * x * std::sqrt(x));
    };

    CentralRectSolver rectSolver;
    TrapezoidSolver trapSolver;
    SimpsonSolver simpSolver;
    GaussSolver gaussSolver;

    Plotter plotter;

    IntegrateTask task_rect(rectSolver);
    IntegrateTask task_trap(trapSolver);
    IntegrateTask task_simp(simpSolver, &plotter);
    IntegrateTask task_gauss(gaussSolver);

    task_rect.run(integrand, a, b, tol);
    task_trap.run(integrand, a, b, tol);
    task_simp.run(integrand, a, b, tol);
    task_gauss.run(integrand, a, b, tol);
}
