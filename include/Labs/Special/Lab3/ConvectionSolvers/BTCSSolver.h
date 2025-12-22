#ifndef NUMERICAL_METHODS_IN_PHYSICS_BTCSSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_BTCSSOLVER_H
#pragma once

#include "Labs/Special/Lab3/Base/IConvectionSolver.h"
#include "Base/ISolver.h"

#include <Eigen/Dense>
#include <vector>
#include <string>

/// Неявная двухслойная схема BTCS для 1D переноса:
/// (U_i^{n+1} - U_i^n)/dt + u * (U_{i+1}^{n+1} - U_{i-1}^{n+1})/(2 dx) = 0
class BTCSSolver : public IConvectionSolver {
public:
    explicit BTCSSolver(ISolver& linearSolver)
        : linSolver(linearSolver) {}

    [[nodiscard]] std::string name() const override { return "BTCS implicit"; }

    void step(std::vector<double>& u_next,
              const std::vector<double>& u,
              double c) const override
    {
        const std::size_t N = u.size();
        if (N < 3) {          // слишком мало точек для внутренней схемы
            u_next = u;
            return;
        }
        if (u_next.size() != N) u_next.resize(N);

        // Количество внутренних узлов i = 1..N-2.
        const std::size_t M = N - 2;

        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(M, M);
        Eigen::VectorXd b(M);

        // Коэффициенты трёхдиагональной BTCS‑матрицы:
        // - (c/2) * U_{i-1}^{n+1} + 1 * U_i^{n+1} + (c/2) * U_{i+1}^{n+1} = U_i^n.
        const double alpha = -0.5 * c;  // поддиагональ
        const double beta  =  1.0;      // диагональ
        const double gamma =  0.5 * c;  // наддиагональ

        for (std::size_t k = 0; k < M; ++k) {
            std::size_t i = k + 1; // глобальный индекс

            // диагональ
            A(k, k) = beta;

            // поддиагональный элемент (i-1 -> k-1)
            if (k > 0)
                A(k, k-1) = alpha;

            // наддиагональный элемент (i+1 -> k+1)
            if (k + 1 < M)
                A(k, k+1) = gamma;

            // правая часть без учёта границ: U_i^n
            b(k) = u[i];
        }

        // Учёт граничных условий (Дирихле: U_0^{n+1} и U_{N-1}^{n+1} считаем известными)
        // Для i = 1 (k = 0): есть член alpha * U_0^{n+1}
        b(0)    -= alpha * u[0];

        // Для i = N-2 (k = M-1): есть член gamma * U_{N-1}^{n+1}
        b(M-1)  -= gamma * u[N-1];

        SolveResult res = linSolver.solve(A, b);
        const Eigen::VectorXd& x = res.solution;

        // Собираем полный вектор u_next
        u_next[0]   = u[0];      // левая граница
        u_next[N-1] = u[N-1];    // правая граница

        for (std::size_t k = 0; k < M; ++k) {
            std::size_t i = k + 1;
            u_next[i] = x(k);
        }
    }

private:
    ISolver& linSolver;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_BTCSSOLVER_H