#ifndef NUMERICAL_METHODS_IN_PHYSICS_MOMENTUMUASSEMBLER_H
#define NUMERICAL_METHODS_IN_PHYSICS_MOMENTUMUASSEMBLER_H
#pragma once
#include <algorithm>
#include <cmath>
#include "Labs/Special/Lab8/Base/StaggeredGrid.h"
#include "Labs/Special/Lab8/Base/Config/SimConfig.h"
#include "Labs/Special/Lab8/Base/Field/ScalarField.h"
#include "Labs/Special/Lab8/Base/Field/UField.h"
#include "Labs/Special/Lab8/Base/Field/VField.h"
#include "Labs/Special/Lab8/Linag/LinearSystem.h"

struct MomentumUAssemblyResult {
    LinearSystem system;
    ScalarField d;
};

class MomentumUAssembler {
public:
    static int unknownCount(const StaggeredGrid& g) {
        return (g.Nx - 1) * g.Ny;
    }

    static int index(const StaggeredGrid& g, int i, int j) {
        return j * (g.Nx - 1) + (i - 1);
    }

    static MomentumUAssemblyResult assemble(const StaggeredGrid& g,
                                            const SimulationConfig& cfg,
                                            const UField& u,
                                            const VField& v,
                                            const ScalarField& p) {
        MomentumUAssemblyResult result;
        result.system.reset(unknownCount(g));
        result.d.resize(g.uNx(), g.uNy(), 0.0);

        const double dx = g.dx();
        const double dy = g.dy();
        const double rho = cfg.rho;
        const double mu = cfg.mu();
        const double eps = 1e-14;

        const double Dw = mu * dy / dx;
        const double De = mu * dy / dx;
        const double Ds = mu * dx / dy;
        const double Dn = mu * dx / dy;

        for (int j = 0; j < g.uNy(); ++j) {
            for (int i = 1; i < g.uNx() - 1; ++i) {
                const int row = index(g, i, j);

                const double uw = 0.5 * (u(i - 1, j) + u(i, j));
                const double ue = 0.5 * (u(i, j) + u(i + 1, j));
                const double vs = 0.5 * (v(i - 1, j) + v(i, j));
                const double vn = 0.5 * (v(i - 1, j + 1) + v(i, j + 1));

                const double Fw = rho * uw * dy;
                const double Fe = rho * ue * dy;
                const double Fs = rho * vs * dx;
                const double Fn = rho * vn * dx;

                double aW = Dw + std::max(Fw, 0.0);
                double aE = De + std::max(-Fe, 0.0);
                double aS = Ds + std::max(Fs, 0.0);
                double aN = Dn + std::max(-Fn, 0.0);

                double rhs = (p(i - 1, j) - p(i, j)) * dy;

                if (i == 1) {
                    rhs += aW * u(0, j);
                    aW = 0.0;
                }
                if (i == g.uNx() - 2) {
                    aE = 0.0;
                }
                if (j == 0) {
                    aS = 0.0;
                }
                if (j == g.uNy() - 1) {
                    aN = 0.0;
                }

                const double wallSouth = (j == 0) ? 2.0 * mu * dx / dy : 0.0;
                const double wallNorth = (j == g.uNy() - 1) ? 2.0 * mu * dx / dy : 0.0;

                const double aP0 = aW + aE + aS + aN
                                 + (Fe - Fw) + (Fn - Fs)
                                 + wallSouth + wallNorth;
                const double aP = aP0 / cfg.alphaU;
                rhs += (aP - aP0) * u(i, j);

                result.system.add(row, row, aP);
                result.system.add_rhs(row, rhs);

                if (i > 1) {
                    result.system.add(row, index(g, i - 1, j), -aW);
                }
                if (i < g.uNx() - 2) {
                    result.system.add(row, index(g, i + 1, j), -aE);
                }
                if (j > 0) {
                    result.system.add(row, index(g, i, j - 1), -aS);
                }
                if (j < g.uNy() - 1) {
                    result.system.add(row, index(g, i, j + 1), -aN);
                }

                const double sumNb = aW + aE + aS + aN;
                result.d(i, j) = dy / std::max(aP - sumNb, eps);
            }
        }

        return result;
    }
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_MOMENTUMUASSEMBLER_H