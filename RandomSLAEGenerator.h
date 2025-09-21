#ifndef NUMERICAL_METHODS_IN_PHYSICS_RANDOMSLAEGENERATOR_H
#define NUMERICAL_METHODS_IN_PHYSICS_RANDOMSLAEGENERATOR_H

#pragma once
#include "IDataGenerator.h"
#include <random>
class RandomSLAEGenerator : public IDataGenerator {

    bool isDiagonallyDominant(const Eigen::MatrixXd& A) {
        for (int i = 0; i < A.rows(); ++i) {
            double diag = std::abs(A(i, i));
            double sum = 0.0;
            for (int j = 0; j < A.cols(); ++j) {
                if (j != i) sum += std::abs(A(i, j));
            }
            if (diag < sum) return false;
        }
        return true;
    }

    Eigen::MatrixXd generateTridiagonalMatrix(size_t n, int low = 1, int high = 10, int diag_high = 30) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> off_diag_dist(low, high);      // для побочных диагоналей
        std::uniform_int_distribution<> diag_dist(diag_high/2, diag_high); // для главной диагонали

        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
        for (size_t i = 0; i < n; ++i) {
            if (i > 0) {
                A(i, i - 1) = off_diag_dist(gen);
            }
            A(i, i) = diag_dist(gen); // главная диагональ специально усиливается
            if (i < n - 1) {
                A(i, i + 1) = off_diag_dist(gen);
            }
        }
        return A;
    }

    // Точное решение x = [0, 1, ..., n]
    Eigen::VectorXd exactSolution(size_t n) const {
        return Eigen::VectorXd::LinSpaced(n, 0, n-1);
    }

public:
    Eigen::MatrixXd generateMatrix(size_t n) override {
        int diagHigh = 30;
        int maxTries = 100;
        for (int attempt = 0; attempt < maxTries; ++attempt) {
            Eigen::MatrixXd matrix = generateTridiagonalMatrix(n, 1, 10, diagHigh);
            if (isDiagonallyDominant(matrix)) {
                return matrix;
            }
            diagHigh += 10; // наращиваем преобладание только при неудаче
        }
        // Если не удалось, явно формируем матрицу с преобладанием
        Eigen::MatrixXd matrix = generateTridiagonalMatrix(n, 1, 10, diagHigh * 2);
        for (size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            if (i > 0)   sum += std::abs(matrix(i, i - 1));
            if (i < n-1) sum += std::abs(matrix(i, i + 1));
            // делаем диагональный элемент строго больше суммы соседних
            matrix(i, i) = sum + 1;
        }
        return matrix;
    }

    // Правая часть b = Ax
    Eigen::VectorXd generateVector(size_t n) override {
        auto A = generateMatrix(n);
        auto x = exactSolution(n);
        return A * x;
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_RANDOMSLAEGENERATOR_H
