#ifndef NUMERICAL_METHODS_IN_PHYSICS_JACOBISOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_JACOBISOLVER_H

#pragma once
#include "Base/ISolver.h"
#include <Eigen/Dense>
#include <vector>

class JacobiSolver : public ISolver {
public:
    using ISolver::ISolver;
    SolveResult solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) override {
        size_t n = A.rows();
        Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd x_new = x;
        std::vector<std::pair<size_t, double>> residuals;
        double b_norm = b.norm();
        double rel_res;
        size_t k = 0;
        for (; k < maxIter.value_or(1000); ++k) {
            for (size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                for (size_t j = 0; j < n; ++j) {
                    if (j != i) sum += A(i, j) * x(j);
                }
                x_new(i) = (b(i) - sum) / A(i, i);
            }
            rel_res = (A * x_new - b).norm() / b_norm;
            residuals.emplace_back(k + 1, rel_res);
            if (rel_res < tolerance.value_or(1e-8))
                return {x_new,k,rel_res, residuals};
            x = x_new;
        }
        return {x_new,k,rel_res, residuals};
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_JACOBISOLVER_H
