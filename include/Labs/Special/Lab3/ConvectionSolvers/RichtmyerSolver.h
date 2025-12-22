#ifndef NUMERICAL_METHODS_IN_PHYSICS_RICHTMYERSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_RICHTMYERSOLVER_H

#pragma once

#include "Labs/Special/Lab3/Base/IConvectionSolver.h"
#include <vector>

class RichtmyerSolver : public IConvectionSolver {
public:
    [[nodiscard]] std::string name() const override { return "Richtmyer (LW two-step)"; }

    void step(std::vector<double>& u_next,
              const std::vector<double>& u,
              double c) const override
    {
        const size_t N = u.size();
        if (u_next.size() != N) u_next.resize(N);
        if (N < 3) { u_next = u; return; }

        std::vector<double> u_star(N);

        // Lax-предиктор: u_i^* = 0.5(U_{i+1}+U_{i-1}) - 0.5*c*(U_{i+1}-U_{i-1})
        for (size_t i = 1; i + 1 < N; ++i)
            u_star[i] = 0.5*(u[i+1] + u[i-1]) - 0.5*c*(u[i+1] - u[i-1]);

        u_star[0]   = u[0];
        u_star[N-1] = u[N-1];

        // U_i^{n+1} = U_i^n - c*(u_{i+1}^* - u_{i-1}^*)
        for (size_t i = 1; i + 1 < N; ++i)
            u_next[i] = u[i] - c * (u_star[i+1] - u_star[i-1]);

        u_next[0]   = u[0];
        u_next[N-1] = u[N-1];
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_RICHTMYERSOLVER_H
