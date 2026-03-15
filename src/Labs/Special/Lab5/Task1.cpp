#include "Helpers/Plotter.h"
#include "Labs/Special/Lab5/Tasks/CylinderHeatTask.h"

#include <iostream>

int main() {
    try {
        Plotter plotter;

        CylinderHeatProblem problem;
        problem.H = 2.0;       // l / a
        problem.q = 5;       // dT/dz|_{z=0} = -q
        problem.theta = 10.0;   // T|_{z=H} = theta
        problem.rhs = [](double r, double z) {
            return 10;
        };

        EigenCGSolver solver(1e-12, 10000);
        CylinderHeatTask task(problem, solver, plotter);

        task.solveSingle(80, 120);

    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }
}