#ifndef NUMERICAL_METHODS_IN_PHYSICS_MOMENTUMVASSEMBLER_H
#define NUMERICAL_METHODS_IN_PHYSICS_MOMENTUMVASSEMBLER_H
#pragma once
#include <algorithm>
#include <cmath>
#include "Labs/Special/Lab8/Base/StaggeredGrid.h"
#include "Labs/Special/Lab8/Base/Config/SimConfig.h"
#include "Labs/Special/Lab8/Base/Field/ScalarField.h"
#include "Labs/Special/Lab8/Base/Field/UField.h"
#include "Labs/Special/Lab8/Base/Field/VField.h"
#include "Labs/Special/Lab8/Linag/LinearSystem.h"

struct MomentumVAssemblyResult {
    LinearSystem system;
    ScalarField d;
};

class MomentumVAssembler {
public:
    static int unknownCount(const StaggeredGrid& g) {
        return g.Nx * (g.Ny - 1);
    }

    static int index(const StaggeredGrid& g, int i, int j) {
        return (j - 1) * g.Nx + i;
    }

    static MomentumVAssemblyResult assemble(const StaggeredGrid& g,
                                            const SimulationConfig& cfg,
                                            const UField& u,
                                            const VField& v,
                                            const ScalarField& p) {
        MomentumVAssemblyResult result;
        result.system.reset(unknownCount(g));
        result.d.resize(g.vNx(), g.vNy(), 0.0);

        const double dx = g.dx();
        const double dy = g.dy();
        const double rho = cfg.rho;
        const double mu = cfg.mu();
        const double eps = 1e-14;

        const double Dw = mu * dy / dx;
        const double De = mu * dy / dx;
        const double Ds = mu * dx / dy;
        const double Dn = mu * dx / dy;

        for (int j = 1; j < g.vNy() - 1; ++j) {
            for (int i = 0; i < g.vNx(); ++i) {
                const int row = index(g, i, j);

                const double uw = 0.5 * (u(i, j - 1) + u(i, j));
                const double ue = 0.5 * (u(i + 1, j - 1) + u(i + 1, j));
                const double vs = 0.5 * (v(i, j - 1) + v(i, j));
                const double vn = 0.5 * (v(i, j) + v(i, j + 1));

                const double Fw = rho * uw * dy;
                const double Fe = rho * ue * dy;
                const double Fs = rho * vs * dx;
                const double Fn = rho * vn * dx;

                double aW = Dw + std::max(Fw, 0.0);
                double aE = De + std::max(-Fe, 0.0);
                double aS = Ds + std::max(Fs, 0.0);
                double aN = Dn + std::max(-Fn, 0.0);

                double rhs = (p(i, j - 1) - p(i, j)) * dx;

                if (i == 0) {
                    rhs += aW * cfg.vInlet;
                    aW = 0.0;
                }
                if (i == g.vNx() - 1) {
                    aE = 0.0;
                }
                if (j == 1) {
                    aS = 0.0;
                }
                if (j == g.vNy() - 2) {
                    aN = 0.0;
                }

                const double aP0 = aW + aE + aS + aN + (Fe - Fw) + (Fn - Fs);
                const double aP = aP0 / cfg.alphaV;
                rhs += (aP - aP0) * v(i, j);

                result.system.add(row, row, aP);
                result.system.add_rhs(row, rhs);

                if (i > 0) {
                    result.system.add(row, index(g, i - 1, j), -aW);
                }
                if (i < g.vNx() - 1) {
                    result.system.add(row, index(g, i + 1, j), -aE);
                }
                if (j > 1) {
                    result.system.add(row, index(g, i, j - 1), -aS);
                }
                if (j < g.vNy() - 2) {
                    result.system.add(row, index(g, i, j + 1), -aN);
                }

                const double sumNb = aW + aE + aS + aN;
                result.d(i, j) = dx / std::max(aP - sumNb, eps);
            }
        }

        return result;
    }
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_MOMENTUMVASSEMBLER_H