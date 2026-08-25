#ifndef NUMERICAL_METHODS_IN_PHYSICS_LINEARSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_LINEARSOLVER_H
#pragma once
#include <vector>
#include "LinearSystem.h"

class LinearSolver {
public:
    virtual ~LinearSolver() = default;
    virtual std::vector<double> solve(const LinearSystem& system) = 0;
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_LINEARSOLVER_H