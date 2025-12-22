#ifndef NUMERICAL_METHODS_IN_PHYSICS_MACCORMACKSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_MACCORMACKSOLVER_H
#pragma once

#include "Labs/Special/Lab3/Base/IConvectionSolver.h"
#include <vector>

class MacCormackSolver : public IConvectionSolver {
public:
    [[nodiscard]] std::string name() const override { return "MacCormack"; }

    void step(std::vector<double>& u_next,
              const std::vector<double>& u,
              double c) const override
    {
        const size_t N = u.size();
        if (u_next.size() != N) u_next.resize(N);
        if (N < 3) { u_next = u; return; }

        std::vector<double> u_pred(N);

        // Предиктор: U_i^p = U_i^n - c*(U_{i+1}^n - U_i^n)
        for (size_t i = 0; i + 1 < N; ++i)
            u_pred[i] = u[i] - c * (u[i+1] - u[i]);
        u_pred[N-1] = u[N-1];

        // Корректор: U_i^{n+1} = 0.5*(U_i^n + U_i^p - c*(U_i^p - U_{i-1}^p))
        for (size_t i = 1; i < N; ++i)
            u_next[i] = 0.5 * (u[i] + u_pred[i]
                               - c * (u_pred[i] - u_pred[i-1]));

        u_next[0] = u[0];
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_MACCORMACKSOLVER_H