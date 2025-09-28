#include "../include/SLAE_ExactMethod_Task.h"
#include "../include/RandomSLAEGenerator.h"
#include "../include/ThomasSolver.h"
#include "../include/AsymptoticFuncs.h"
#include "../include/SORSolver.h"
#include <vector>

int main() {
    Plotter plot;
    RandomSLAEGenerator gen;
    /*ThomasSolver solver;

    SLAE_ExactMethod_Task lab(gen, solver, plot, SLAE_ExactMethod_Task::ExpectedCurveFunc(linearAsymptotic));

    std::vector<size_t> sizes;
    for (size_t n = 1000; n <= 4000; n += 100)
        sizes.push_back(n);
    lab.run(sizes);*/

    std::vector<double> omega_values;
    std::vector<double> iteration_counts;
    size_t n = 100;
    size_t min_iters = std::numeric_limits<size_t>::max();
    double omega_opt = 1.0;
    auto* slaeGen = dynamic_cast<RandomSLAEGenerator*>(&gen);
    Eigen::MatrixXd A = slaeGen->generateMatrix(n, true);
    Eigen::VectorXd x_exact = slaeGen->exactSolution(n);
    Eigen::VectorXd b = A * x_exact;
    SORSolver solver(1000,1e-8);

    for (double omega = 0.01; omega < 2.0; omega += 0.01) {

        solver.param = omega;

        auto result = solver.solve(A, b);
        // Найти номер итерации, где невязка стала меньше tol
        size_t iter_to_converge = result.residuals.size();
        omega_values.push_back(omega);
        iteration_counts.push_back(iter_to_converge);

        if (iter_to_converge < min_iters) {
            min_iters = iter_to_converge;
            omega_opt = omega;
        }
    }
    plot.plot(omega_values, iteration_counts, "SOR: iterations vs omega", "omega", "iterations", true);
    std::cout<<"Optimal omega is: "<<omega_opt;



}
