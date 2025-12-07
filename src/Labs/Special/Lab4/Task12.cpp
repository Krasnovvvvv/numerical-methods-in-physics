#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#include "../include/Labs/Lab3/IntegralSolvers/GaussSolver.h"
#include "../include/Labs/Lab3/IntegralSolvers/TrapezoidSolver.h"
#include "../include/Labs/Special/Lab4/Tasks/IntegralEquationTask.h"
#include "Helpers/Plotter.h"

// Fredholm: x(t) + ∫_1^2 (t + 2s) x(s) ds = t^2 - 1
// точное решение: t^2 - 5/32 t - 373/192
void run_fredholm_task(Plotter& plotter, unsigned short mode) {
    const double a = 1.0;
    const double b = 2.0;
    const std::size_t n = 10;
    const double lambda = 1.0;

    TrapezoidSolver integrator;
    IntegralEquationTask task(integrator, &plotter, mode);

    IntegralEquationTask::Kernel kernel = [](double t, double s) {
        return t + 2.0 * s;
    };
    IntegralEquationTask::RHS rhs = [](double t) {
        return t * t - 1.0;
    };
    IntegralEquationTask::Solution exact = [](double t) {
        return t * t - (5.0 / 32.0) * t - 373.0 / 192.0;
    };

    auto x_num = task.solve_fredholm(lambda, a, b, n, kernel, rhs, exact);

    std::cout << "\nFredholm solution at grid nodes:\n";
    double h = (b - a) / static_cast<double>(n - 1);
    for (std::size_t i = 0; i < x_num.size(); ++i) {
        double t  = a + h * static_cast<double>(i);
        double xe = exact(t);
        std::cout << "t = " << t
                  << "  x_num = " << std::fixed << std::setprecision(8) << x_num[i]
                  << "  x_exact = " << xe
                  << "  error = " << std::abs(x_num[i] - xe)
                  << '\n';
    }
}

// Volterra: x(t) = 4 ∫_0^t (s - t) x(s) ds + 3 sin t,
// t ∈ [0,5], τ = 0.25; точное решение: 2 sin(2t) - sin t
void run_volterra_task(Plotter& plotter, unsigned short mode) {
    const double a = 0.0;
    const double b = 5.0;
    const double tau = 0.25;
    const std::size_t n =
        static_cast<std::size_t>((b - a) / tau) + 1;

    TrapezoidSolver integrator;
    IntegralEquationTask task(integrator, &plotter, mode);

    IntegralEquationTask::Kernel kernel = [](double t, double s) {
        return (s - t);
    };
    IntegralEquationTask::RHS rhs = [](double t) {
        return 3.0 * std::sin(t);
    };
    IntegralEquationTask::Solution exact = [](double t) {
        return 2.0 * std::sin(2.0 * t) - std::sin(t);
    };

    auto x_num = task.solve_volterra(a, b, n, kernel, rhs, exact, 4.0);

    std::cout << "\nVolterra solution on grid with step 0.25:\n";
    for (std::size_t i = 0; i < x_num.size(); ++i) {
        double t  = a + tau * static_cast<double>(i);
        double xe = exact(t);
        std::cout << "t = " << t
                  << "  x_num = " << std::fixed << std::setprecision(8) << x_num[i]
                  << "  x_exact = " << xe
                  << "  error = " << std::abs(x_num[i] - xe)
                  << '\n';
    }
}

int main() {
    Plotter plotter;

    unsigned short mode = 1;
    run_fredholm_task(plotter, mode);
    run_volterra_task(plotter, mode);

}
