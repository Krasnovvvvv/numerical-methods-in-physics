#include "Labs/Special/Lab5/Tasks/FEMCylinderHeatTask.h"

int main() {
    try {
        CylinderHeatProblem problem;
        problem.H = 2.0;       // l / a
        problem.q = 5;       // dT/dz|_{z=0} = -q
        problem.theta = 10.0;   // T|_{z=H} = theta
        problem.rhs = [](double r, double z) {
            return 10;
        };

        EigenCGSolver solver(1e-12, 10000);

        FEMCylinderHeatTask task(problem, solver);
        task.solveSingle(80, 120);

    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }
}