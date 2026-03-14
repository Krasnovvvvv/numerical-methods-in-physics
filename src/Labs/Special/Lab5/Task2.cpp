#include "Labs/Special/Lab5/Tasks/FEMCylinderHeatTask.h"

int main() {
    try {
        CylinderHeatProblem problem;
        problem.H = 1.0;       // l / a
        problem.q = 65;       // dT/dz|_{z=0} = -q
        problem.theta = 5.0;   // T|_{z=H} = theta
        problem.rhs = [](double r, double z) {
            return 0;
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