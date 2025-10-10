#include "Base/IIntegralSolver.h"
#include "Labs/Lab3/IntegralSolvers/CentralRectSolver.h"
#include "Labs/Lab3/IntegralSolvers/TrapezoidSolver.h"
#include "Labs/Lab3/IntegralSolvers/SimpsonSolver.h"
#include "Labs/Lab3/Tasks/IntegrateTask.h"
#include <cmath>

int main() {
    constexpr double a = 0.1, b = 0.4, tol = 1e-6;

    auto integrand = [](double x) {
        return std::atan(0.3 * std::pow(x, 4) - 4.0 * x * std::sqrt(x));
    };

    CentralRectSolver rect;
    TrapezoidSolver trap;
    SimpsonSolver simp;

    IntegrateTask task_rect(rect, "The middle rectangles");
    IntegrateTask task_trap(trap, "Trapezoids");
    IntegrateTask task_simp(simp, "Simpson");

    task_rect.run(integrand, a, b, tol);
    task_trap.run(integrand, a, b, tol);
    task_simp.run(integrand, a, b, tol);
}
