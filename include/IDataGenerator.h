
#ifndef NUMERICAL_METHODS_IN_PHYSICS_IDATAGENERATOR_H
#define NUMERICAL_METHODS_IN_PHYSICS_IDATAGENERATOR_H

#pragma once
#include <Eigen/Dense>
class IDataGenerator {
public:
    virtual Eigen::MatrixXd generateMatrix(size_t n, bool is_tridiagonal = false) = 0;
    virtual Eigen::VectorXd generateVector(size_t n) = 0;
    virtual ~IDataGenerator() = default;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_IDATAGENERATOR_H
