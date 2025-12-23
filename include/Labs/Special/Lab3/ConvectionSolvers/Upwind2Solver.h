#ifndef NUMERICAL_METHODS_IN_PHYSICS_UPWIND2SOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_UPWIND2SOLVER_H

#pragma once

#include "Labs/Special/Lab3/Base/IConvectionSolver.h"

#include <vector>
#include <string>
#include <cstddef>

class Upwind2Solver : public IConvectionSolver {
public:
    [[nodiscard]] std::string name() const override {
        return "Upwind 2nd order (with 1st-order boundary)";
    }

    void step(std::vector<double>&       u_next,
              const std::vector<double>& u,
              double                     c) const override
    {
        const std::size_t N = u.size();
        if (u_next.size() != N)
            u_next.resize(N);

        if (N < 3) {
            u_next = u;
            return;
        }

        if (c >= 0.0) {
            // -------- u > 0: волна вправо, левый upwind --------

            // 1) Первый внутренний узел: upwind 1-го порядка (i = 1)
            if (N > 2) {
                double upwind1 = u[1] - u[0];
                u_next[1] = u[1] - c * upwind1;
            }

            // 2) Узлы второго порядка: i = 2..N-2
            for (std::size_t i = 2; i + 1 < N; ++i) {
                double upwind1 = u[i] - u[i - 1];
                double corr    = u[i] - 2.0 * u[i - 1] + u[i - 2];

                u_next[i] = u[i]
                            - c * upwind1
                            - 0.5 * c * (1.0 - c) * corr;
            }

            // Узлы 0 и N-1 остаются под управлением ConvectionTask.
        } else {
            // -------- u < 0: волна влево, правый upwind --------

            // 1) Первый внутренний узел слева от правой границы: i = N-2
            //    upwind 1-го порядка: U_{N-2}^{n+1} = U_{N-2}^n - c (U_{N-2}^n - U_{N-1}^n)
            if (N > 2) {
                std::size_t i = N - 2;
                double upwind1 = u[i] - u[i + 1];
                u_next[i] = u[i] - c * upwind1;
            }

            // 2) Узлы второго порядка: i = 1..N-3
            for (std::size_t i = 1; i + 2 < N; ++i) {
                double upwind1 = u[i] - u[i + 1];
                double corr    = u[i] - 2.0 * u[i + 1] + u[i + 2];

                u_next[i] = u[i]
                            - c * upwind1
                            - 0.5 * c * (1.0 - c) * corr;
            }

            // Узлы 0 и N-1 снова задаются только в ConvectionTask.
        }
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_UPWIND2SOLVER_H
