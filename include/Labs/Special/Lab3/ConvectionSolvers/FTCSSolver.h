#ifndef NUMERICAL_METHODS_IN_PHYSICS_FTCSSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_FTCSSOLVER_H

#pragma once

#include "Labs/Special/Lab3/Base/IConvectionSolver.h"

class FTCSSolver : public IConvectionSolver {
public:
    [[nodiscard]] std::string name() const override { return "FTCS"; }

    void step(std::vector<double>& u_next,
              const std::vector<double>& u,
              double c) const override
    {
        const size_t N = u.size();
        if (u_next.size() != N) u_next.resize(N);

        // U_i^{n+1} = U_i^n - (c/2)*(U_{i+1}^n - U_{i-1}^n)
        for (size_t i = 1; i + 1 < N; ++i)
            u_next[i] = u[i] - 0.5 * c * (u[i+1] - u[i-1]);

        if (N >= 1) u_next[0]   = u[0];
        if (N >= 2) u_next[N-1] = u[N-1];
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_FTCSSOLVER_H
