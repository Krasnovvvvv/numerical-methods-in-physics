#ifndef NUMERICAL_METHODS_IN_PHYSICS_UPWIND2SOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_UPWIND2SOLVER_H

#pragma once

#include "Labs/Special/Lab3/Base/IConvectionSolver.h"

#include <vector>
#include <string>
#include <cmath>

class Upwind2Solver : public IConvectionSolver {

public:
    [[nodiscard]] std::string name() const override {
        return "Upwind 2nd order";
    }

    void step(std::vector<double>& u_next,
              const std::vector<double>& u,
              double c) const override
    {
        const size_t N = u.size();
        if (u_next.size() != N) u_next.resize(N);
        if (N < 3) { u_next = u; return; }

        // u>0: используем узлы i, i-1, i-2
        if (c >= 0.0) {

            u_next[0] = u[0];
            u_next[1] = u[1];

            for (size_t i = 2; i < N; ++i) {
                double upwind1 = u[i] - u[i-1];
                double corr    = u[i] - 2.0*u[i-1] + u[i-2];

                u_next[i] = u[i]
                            - c * upwind1
                            - 0.5 * c * (1.0 - c) * corr;
            }
        }
        // u<0: зеркальная формула, используем i, i+1, i+2
        else {
            u_next[N-1] = u[N-1];
            u_next[N-2] = u[N-2];

            for (size_t i = 0; i + 2 < N; ++i) {
                double upwind1 = u[i] - u[i+1];
                double corr    = u[i] - 2.0*u[i+1] + u[i+2];

                u_next[i] = u[i]
                            - c * upwind1
                            - 0.5 * c * (1.0 - c) * corr;
            }
        }
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_UPWIND2SOLVER_H
