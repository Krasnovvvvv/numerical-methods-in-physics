#include "Labs/Lab1/SLAEGenerators/RandomSLAEGenerator.h"
#include "Labs/Lab1/SLAESolvers/SORSolver.h"
#include "Helpers/Plotter.h"
#include <vector>

int main() {

    Plotter plot;
    RandomSLAEGenerator gen;
    auto* slaeGen = dynamic_cast<RandomSLAEGenerator*>(&gen);
    std::vector<double> omega_values;
    std::vector<double> iteration_counts;

    size_t n = 100;
    double omega_opt = 1.0;
    size_t min_iters = std::numeric_limits<size_t>::max();

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
    plot.plot(omega_values, iteration_counts,
              "SOR: iterations vs omega",
              "omega", "iterations", true);
    std::cout<<"Optimal omega is: "<<omega_opt;



}

