#include "Labs/Lab1/SLAESolvers/ThomasSolver.h"
#include "Labs/Lab5/Tasks/TaskBessel.h"
#include <cmath>

int main () {

    size_t N = 100;
    double a = 1.0, b = 10.0;
    int order = 1;
    double mu_left = 1.0, mu_right = 0.0;

    ThomasSolver solver;
    Plotter plotter;

    TaskBessel task(solver, &plotter);
    task.run(order, a, b, mu_left, mu_right, N);
}