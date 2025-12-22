#ifndef NUMERICAL_METHODS_IN_PHYSICS_LAXWENDROFFSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_LAXWENDROFFSOLVER_H

#pragma once

#include "Labs/Special/Lab3/Base/IConvectionSolver.h"

class LaxWendroffSolver : public IConvectionSolver {
public:
    [[nodiscard]] std::string name() const override { return "Lax–Wendroff one-step"; }

    void step(std::vector<double>& u_next,
              const std::vector<double>& u,
              double c) const override
    {
        const size_t N = u.size();
        if (u_next.size() != N) u_next.resize(N);
        const double c2 = c * c;

        // U_i^{n+1} = U_i^n - (c/2)(U_{i+1}-U_{i-1})
        //             + (c^2/2)(U_{i+1}-2U_i+U_{i-1})
        for (size_t i = 1; i + 1 < N; ++i) {
            u_next[i] = u[i]
                      - 0.5 * c  * (u[i+1] - u[i-1])
                      + 0.5 * c2 * (u[i+1] - 2.0*u[i] + u[i-1]);
        }

        if (N >= 1) u_next[0]   = u[0];
        if (N >= 2) u_next[N-1] = u[N-1];
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_LAXWENDROFFSOLVER_H
