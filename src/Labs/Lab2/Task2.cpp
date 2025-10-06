#include "Labs/Lab2/RootSolvers/BruteForceBisecRootSolver.h"
#include "Labs/Lab2/RootSolvers/ParallelBisecRootSolver.h"
#include "Labs/Lab2/RootSolvers/SimpleIterSolver.h"
#include "Labs/Lab2/RootSolvers/NewtonSolver.h"
#include "Labs/Lab2/Tasks/LampVACHTask.h"
#include "Helpers/Plotter.h"

int main() {
    double R0 = 150.0;
    double alpha = 4.5e-3;
    double T0 = 293.0;
    double sigma = 5.67e-8;
    double S = 0.001;
    double epsilon = 0.3;
    double V_min = 0.1, V_max = 220.0, dV = 0.2;
    double tol = 1e-10;
    size_t max_iter = 100;

    Plotter plotter;

    BruteForceBisecRootSolver brute_solver;
    LampVACHTask brute_task(brute_solver, plotter);
    brute_task.run(V_min, V_max, dV, R0, alpha, T0, sigma, S, epsilon, tol, max_iter, 500, "Brute Force");

    ParallelBisecRootSolver parallel_solver(6);
    LampVACHTask parallel_task(parallel_solver, plotter);
    parallel_task.run(V_min, V_max, dV, R0, alpha, T0, sigma, S, epsilon, tol, max_iter, 500, "Parallel");

    NewtonSolver newton_solver;
    LampVACHTask newton_task(newton_solver, plotter);
    newton_task.run(V_min, V_max, dV, R0, alpha, T0, sigma, S, epsilon, tol, max_iter, 500, "Newton");

    SimpleIterSolver iter_solver(0.25);
    LampVACHTask iter_task(iter_solver, plotter);
    iter_task.run(V_min, V_max, dV, R0, alpha, T0, sigma, S, epsilon, tol, max_iter, 500, "Simple Iteration");


}
