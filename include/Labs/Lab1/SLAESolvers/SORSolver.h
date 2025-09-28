#ifndef NUMERICAL_METHODS_IN_PHYSICS_SORSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_SORSOLVER_H

#pragma once
#include "Base/ISolver.h"
#include <Eigen/Dense>
#include <vector>

class SORSolver : public ISolver {
public:
    using ISolver::ISolver;

    SolveResult solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) override {
        size_t n = A.rows();
        Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
        double b_norm = b.norm();
        std::vector<std::pair<size_t, double>> residuals;
        double rel_res = 0.0;
        size_t k = 0;
        for (; k < maxIter.value_or(1000); ++k) {
            for (size_t i = 0; i < n; ++i) {
                double sum1 = 0.0, sum2 = 0.0;
                for (size_t j = 0; j < i; ++j) sum1 += A(i, j) * x(j);
                for (size_t j = i + 1; j < n; ++j) sum2 += A(i, j) * x(j);
                double x_new = (b(i) - sum1 - sum2) / A(i, i);
                x(i) = (1 - param.value_or(1.3)) * x(i) + param.value_or(1.3) * x_new;
            }
            rel_res = (A * x - b).norm() / b_norm;
            residuals.emplace_back(k + 1, rel_res);
            if (rel_res < tolerance.value_or(1e-8))
                return {x, k + 1, rel_res, residuals};
        }
        return {x, k, rel_res, residuals};
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_SORSOLVER_H
