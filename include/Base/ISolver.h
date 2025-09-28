#ifndef NUMERICAL_METHODS_IN_PHYSICS_ISOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_ISOLVER_H

#pragma once
#include <Eigen/Dense>
#include <optional>

struct SolveResult {
    Eigen::VectorXd solution;
    size_t iterations;
    double residual;
    std::vector<std::pair<size_t, double>> residuals;
};

class ISolver {
public:
    virtual SolveResult solve(const Eigen::MatrixXd&, const Eigen::VectorXd&) = 0;
    virtual ~ISolver() = default;

    explicit ISolver(std::optional<int> max_iter = std::nullopt,
            std::optional<double> tol = std::nullopt,
            std::optional<double> prm = std::nullopt) : maxIter(max_iter), tolerance(tol), param(prm) {}

    std::optional<double> param;

protected:
    std::optional<int> maxIter;
    std::optional<double> tolerance;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_ISOLVER_H
