#ifndef NUMERICAL_METHODS_IN_PHYSICS_ISOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_ISOLVER_H

#pragma once
#include <Eigen/Dense>
struct SolveResult {
    Eigen::VectorXd solution;
    size_t iterations;
    double residual;
};

class ISolver {
public:
    virtual SolveResult solve(const Eigen::MatrixXd&, const Eigen::VectorXd&) = 0;
    virtual ~ISolver() = default;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_ISOLVER_H
