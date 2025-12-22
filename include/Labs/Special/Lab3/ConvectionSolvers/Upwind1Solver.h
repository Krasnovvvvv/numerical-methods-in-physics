#ifndef NUMERICAL_METHODS_IN_PHYSICS_UPWIND1SOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_UPWIND1SOLVER_H
#pragma once

#include "Labs/Special/Lab3/Base/IConvectionSolver.h"

class Upwind1Solver : public IConvectionSolver {
public:
    [[nodiscard]] std::string name() const override { return "Upwind 1st order"; }

    void step(std::vector<double>& u_next,
              const std::vector<double>& u,
              double c) const override
    {
        const size_t N = u.size();
        if (u_next.size() != N) u_next.resize(N);
        if (N < 2) { u_next = u; return; }

        if (c >= 0.0) {
            // U_i^{n+1} = U_i^n - c*(U_i^n - U_{i-1}^n)
            u_next[0] = u[0];
            for (size_t i = 1; i < N; ++i)
                u_next[i] = u[i] - c * (u[i] - u[i-1]);
        } else {
            // вариант для u<0: U_i^{n+1} = U_i^n - c*(U_{i+1}^n - U_i^n)
            for (size_t i = 0; i + 1 < N; ++i)
                u_next[i] = u[i] - c * (u[i+1] - u[i]);
            u_next[N-1] = u[N-1];
        }
    }
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_UPWIND1SOLVER_H