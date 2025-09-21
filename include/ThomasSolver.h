#ifndef NUMERICAL_METHODS_IN_PHYSICS_THOMASSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_THOMASSOLVER_H

#pragma once
#include "ISolver.h"
#include <Eigen/Dense>
#include <vector>

class ThomasSolver : public ISolver {
public:
    SolveResult solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) override {
        size_t n = b.size();
        std::vector<double> a(n), b_diag(n), c(n), f(n);
        for (size_t i = 0; i < n; ++i) {
            if (i > 0)   a[i] = -A(i, i-1);
            if (i < n-1) b_diag[i] = -A(i, i+1);
            c[i] = A(i, i);
            f[i] = b[i];
        }

        std::vector<double> alpha(n), beta(n);
        alpha[0] = b_diag[0] / c[0];
        beta[0]  = f[0] / c[0];

        for (size_t i = 1; i < n-1; ++i) {
            double denom = c[i] - a[i] * alpha[i-1];
            alpha[i] = b_diag[i] / denom;
            beta[i]  = (f[i] + a[i] * beta[i-1]) / denom;
        }
        beta[n-1] = (f[n-1] + a[n-1] * beta[n-2]) / (c[n-1] - a[n-1] * alpha[n-2]);

        Eigen::VectorXd x(n);
        x[n-1] = beta[n-1];
        for (int i = n-2; i >= 0; --i)
            x[i] = alpha[i] * x[i+1] + beta[i];

        double residual = (A * x - Eigen::Map<const Eigen::VectorXd>(f.data(), n)).norm();
        return {x, 0, residual};
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_THOMASSOLVER_H
