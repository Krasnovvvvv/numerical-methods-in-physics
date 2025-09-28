#ifndef NUMERICAL_METHODS_IN_PHYSICS_CONJUGATEGRADIENTSSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_CONJUGATEGRADIENTSSOLVER_H

#pragma once
#include "ISolver.h"
#include <Eigen/Dense>
#include <vector>

class ConjugateGradientsSolver : public ISolver {
public:
    using ISolver::ISolver;

    SolveResult solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) override {
        size_t n = A.rows();
        Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd r = b - A * x;
        Eigen::VectorXd s = r;
        double b_norm = b.norm();
        std::vector<std::pair<size_t, double>> residuals;
        size_t k = 0;
        for (; k < maxIter.value_or(1000); ++k) {
            double alpha = r.dot(r) / s.dot(A * s);
            x = x + alpha * s;
            Eigen::VectorXd r_new = r - alpha * A * s;
            double rel_res = r_new.norm() / b_norm;
            residuals.emplace_back(k + 1, rel_res);
            if (rel_res < tolerance.value_or(1e-8))
                return {x, k + 1, rel_res, residuals};
            double beta = r_new.dot(r_new) / r.dot(r);
            s = r_new + beta * s;
            r = r_new;
        }
        return {x, k, (A * x - b).norm() / b_norm, residuals};
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_CONJUGATEGRADIENTSSOLVER_H
