#ifndef NUMERICAL_METHODS_IN_PHYSICS_EIGENCGSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_EIGENCGSOLVER_H

#pragma once

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <stdexcept>

#include "Labs/Special/Lab5/Base/SolveResult.h"

class EigenCGSolver {
public:
    explicit EigenCGSolver(double tolerance = 1e-12, int maxIterations = 10000)
        : tolerance_(tolerance), maxIterations_(maxIterations) {}

    [[nodiscard]] SolveResult solve(const Eigen::SparseMatrix<double>& A,
                      const Eigen::VectorXd& b) const
    {
        Eigen::ConjugateGradient<
            Eigen::SparseMatrix<double>,
            Eigen::Lower | Eigen::Upper
        > solver;

        solver.setTolerance(tolerance_);
        solver.setMaxIterations(maxIterations_);
        solver.compute(A);

        if (solver.info() != Eigen::Success)
            throw std::runtime_error("Eigen CG compute() failed.");

        Eigen::VectorXd x = solver.solve(b);

        if (solver.info() != Eigen::Success)
            throw std::runtime_error("Eigen CG solve() failed.");

        SolveResult result;
        result.solution = x;
        result.iterations = solver.iterations();
        result.estimatedError = solver.error();
        result.relativeResidual = (A * x - b).norm() / b.norm();

        return result;
    }

private:
    double tolerance_;
    int maxIterations_;
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_EIGENCGSOLVER_H
