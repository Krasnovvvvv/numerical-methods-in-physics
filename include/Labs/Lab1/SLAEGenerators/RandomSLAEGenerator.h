#ifndef NUMERICAL_METHODS_IN_PHYSICS_RANDOMSLAEGENERATOR_H
#define NUMERICAL_METHODS_IN_PHYSICS_RANDOMSLAEGENERATOR_H

#pragma once
#include "Base/IDataGenerator.h"
#include <random>
#include <climits>

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

public:

    explicit RandomSLAEGenerator(int min_v = 10, int max_v = 100) : minValue(min_v), maxValue(max_v) {
        if (minValue < INT_MIN) minValue = INT_MIN;
        if (minValue > INT_MAX) minValue = INT_MAX;
        if (maxValue < INT_MIN) maxValue = INT_MIN;
        if (maxValue > INT_MAX) maxValue = INT_MAX;
    }

    Eigen::MatrixXd generateMatrix(size_t n, bool is_tridiagonal = false) override {
        if(is_tridiagonal) {
            int diagHigh = 2*maxValue;
            int maxTries = 100;
            for (int attempt = 0; attempt < maxTries; ++attempt) {
                Eigen::MatrixXd matrix = generateTridiagonalMatrix(n, minValue, maxValue, diagHigh);
                if (isDiagonallyDominant(matrix)) {
                    return matrix;
                }
                diagHigh += 10;
            }
            // Если не удалось, явно формируем матрицу с преобладанием
            Eigen::MatrixXd matrix = generateTridiagonalMatrix(n, minValue, maxValue, diagHigh * 2);
            for (size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                if (i > 0) sum += std::abs(matrix(i, i - 1));
                if (i < n - 1) sum += std::abs(matrix(i, i + 1));
                // делаем диагональный элемент строго больше суммы соседних
                matrix(i, i) = sum + 1;
            }
            return matrix;
        }
        else {

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(minValue, maxValue);
            Eigen::MatrixXd matrix(n, n);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    matrix(i, j) = dis(gen);
                }
            }
            return matrix;
        }
    }

    // Точное решение x = [0, 1, ..., n-1]
    Eigen::VectorXd exactSolution(size_t n) const {
        return Eigen::VectorXd::LinSpaced(n, 0, n-1);
    }

    Eigen::VectorXd generateVector(size_t n) override {
        throw std::logic_error("Not used in this context");
    }

    int minValue;
    int maxValue;

};

#endif //NUMERICAL_METHODS_IN_PHYSICS_RANDOMSLAEGENERATOR_H
