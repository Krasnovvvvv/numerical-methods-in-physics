#include "Base/IRootSolver.h"
#include "Labs/Lab2/RootSolvers/BruteForceBisecRootSolver.h"
#include "Labs/Lab2/RootSolvers/ParallelBisecRootSolver.h"
#include "Labs/Lab2/RootSolvers/SimpleIterSolver.h"
#include "Labs/Lab2/RootSolvers/NewtonSolver.h"
#include "Labs/Lab2/Tasks/RootFindTask.h"
#include <cmath>
#include <iostream>

int main() {
    auto func = [](double x) {
        return 2 * std::log(x) + std::sin(std::log(x)) - std::cos(std::log(x));
    };
    auto isInDomain = [](double x) { return x > 0; };

    double a = 1, b = 3;
    double tol = 1e-12;
    double step = 1e-5;
    size_t max_iter = 5000;

    std::cout << "Brute Force Solver:\n";
    BruteForceBisecRootSolver brute_solver;
    RootFindTask brute_task(brute_solver);
    brute_task.run(func, isInDomain, tol, a, b, step, max_iter);

    std::cout << "\nParallel Solver:\n";
    ParallelBisecRootSolver parallel_solver(6);
    RootFindTask parallel_task(parallel_solver);
    parallel_task.run(func, isInDomain, tol, a, b, step, max_iter);

    std::cout << "\nSimple Iteration Solver:\n";
    SimpleIterSolver iter_solver(0.8);
    RootFindTask iter_task(iter_solver);
    iter_task.run(func, isInDomain, tol, 1.4, 0, 0, max_iter);

    std::cout << "\nNewton Solver:\n";
    NewtonSolver newton_solver;
    RootFindTask newton_task(newton_solver);
    newton_task.run(func, isInDomain, tol, 1.4, 0, 0, max_iter);

    return 0;
}
