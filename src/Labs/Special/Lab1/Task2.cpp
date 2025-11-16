#include "Labs/Special/Lab1/Tasks/ODETask.h"
#include "Labs/Lab4/ODESolvers/RungeKutta4Solver.h"
#include "Labs/Special/Lab1/ODESolvers/GearBDFSolver.h"
#include "Labs/Special/Lab1/ODESolvers/AdamsMoulton4Solver.h"
#include "Labs/Special/Lab1/ODEs/ODEs.h"
#include "Helpers/Plotter.h"
#include <cmath>
#include <functional>
#include <iostream>

void run_corrections_sweep(
    std::function<std::vector<double>(double, const std::vector<double>&)> func,
    std::vector<double> y0, double t0, double tn, double h,
    std::function<double(double)> exact_sol,
    Plotter* plotter,
    int min_corr, int max_corr
) {
    RungeKutta4Solver rk4solver;
    std::vector<std::vector<double>> xs, ys;
    std::vector<std::string> labels;

    for (int corr = min_corr; corr <= max_corr; ++corr) {
        AdamsMoulton4Solver solver(&rk4solver, corr);
        ODEResult result = solver.solve(func, y0, t0, tn, h);

        std::vector<double> x, y_err;
        for (size_t i = 0; i < result.t.size(); ++i) {
            x.push_back(result.t[i]);
            double exact = exact_sol ? exact_sol(result.t[i]) : NAN;
            if (exact_sol)
                y_err.push_back(std::abs(result.y[i][0] - exact));
        }
        xs.push_back(x);
        ys.push_back(y_err);
        labels.push_back("corr=" + std::to_string(corr));
        std::cout << "Max error for " << corr << " corrections: "
                  << (y_err.empty() ? 0.0 : *std::max_element(y_err.begin(), y_err.end())) << std::endl;
    }
    // Draw all the errors on one graph
    plotter->plot(xs, ys, labels, "t", "Error");
}

int main() {

    // === The initial task ===
    std::vector<double> y0 = {0.5};
    constexpr double t0 = 0.0, tn = 0.1, h = 0.001;

    auto U_exact = [](double t) {
        double sqrt05 = std::sqrt(0.5);
        double phi0 = std::atan(sqrt05 * 0.5);
        return (1.0 / sqrt05) * std::tan(sqrt05 * t + phi0);
    };

    Plotter plot;
    RungeKutta4Solver rk4solver;

    // --- Gear method of the 2nd order with RK4 start (numerical and accurate graph) ---
    std::cout << "\nGear 2 with RK4 start:\n";
    {
        GearBDFSolver solver1(2, &rk4solver);

        // --- numeric vs analytic ---
        special::ODETask task1(solver1, &plot, 1, {"U"});
        task1.run(ode_rhs, y0, t0, tn, h, {U_exact});

        // --- |y_exact - y_num| ---
        special::ODETask task2(solver1, &plot, 2, {"U"});
       task2.run(ode_rhs, y0, t0, tn, h, {U_exact});
    }

    // --- ABM4: numerical and accurate, 1 correction ---
    std::cout << "\nAdams-Bashforth-Moulton 4 (1 correction):\n";
    {
        AdamsMoulton4Solver abm4solver1(&rk4solver, 1);
        special::ODETask abm_task1(abm4solver1, &plot, 1, {"U"});
        abm_task1.run(ode_rhs, y0, t0, tn, h, {U_exact});

        special::ODETask abm_task1_err(abm4solver1, &plot, 2, {"U"});
        abm_task1_err.run(ode_rhs, y0, t0, tn, h, {U_exact});
    }

    // --- Investigation of the error dependence on the number of corrections ---
    std::cout << "\nSweep by the number of correction iterations in ABM4:\n";
    run_corrections_sweep(ode_rhs, y0, t0, tn, h, U_exact, &plot, 1, 5);

    // --- ADDITIONAL BLOCK: A Tough task ---
    std::cout << "\n=== A tough task for a correction test ===\n";
    {
        auto rhs_stiff = [](double t, const std::vector<double>& y) {
            return std::vector<double>{ -1000.0*y[0] + 3000.0 - 2000.0*std::exp(-t) };
        };
        auto exact_stiff = [](double t) {
            return 3. - 0.998 * std::exp(-1000.0 * t) - 2.002 * std::exp(-t);
        };
        std::vector<double> y0_stiff = {0.0};
        double t0_stiff = 0.0, tn_stiff = 0.01, h_stiff = 0.001;

        std::cout << "Sweep (corr=1..) for stiff problem, h=" << h_stiff << "\n";
        run_corrections_sweep(rhs_stiff, y0_stiff, t0_stiff, tn_stiff, h_stiff,
            exact_stiff, &plot, 1, 10);
    }
}
