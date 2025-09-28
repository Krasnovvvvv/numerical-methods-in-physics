#ifndef NUMERICAL_METHODS_IN_PHYSICS_GRADIENTDESCENTSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_GRADIENTDESCENTSOLVER_H

#pragma once
#include "ISolver.h"
#include <Eigen/Dense>
#include <vector>

class GradientDescentSolver : public ISolver {
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
            Eigen::VectorXd r = b - A * x;
            double alpha = r.dot(r) / r.dot(A * r);
            x = x + alpha * r;
            rel_res = (A * x - b).norm() / b_norm;
            residuals.emplace_back(k + 1, rel_res);
            if (rel_res < tolerance.value_or(1e-8))
                return {x, k + 1, rel_res, residuals};
        }
        return {x, k, rel_res, residuals};
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_GRADIENTDESCENTSOLVER_H
