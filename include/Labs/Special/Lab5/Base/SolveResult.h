#ifndef NUMERICAL_METHODS_IN_PHYSICS_SOLVERESULT_H
#define NUMERICAL_METHODS_IN_PHYSICS_SOLVERESULT_H

#pragma once
#include <Eigen/Dense>

struct SolveResult {
    Eigen::VectorXd solution;
    int iterations = 0;
    double estimatedError = 0.0;
    double relativeResidual = 0.0;
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_SOLVERESULT_H
