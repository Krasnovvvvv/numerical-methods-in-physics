#ifndef NUMERICAL_METHODS_IN_PHYSICS_PRESSURECORRECTIONASSEMBLER_H
#define NUMERICAL_METHODS_IN_PHYSICS_PRESSURECORRECTIONASSEMBLER_H
#pragma once
#include "Labs/Special/Lab8/Base/StaggeredGrid.h"
#include "Labs/Special/Lab8/Base/Config/SimConfig.h"
#include "Labs/Special/Lab8/Base/Field/ScalarField.h"
#include "Labs/Special/Lab8/Base/Field/UField.h"
#include "Labs/Special/Lab8/Base/Field/VField.h"
#include "Labs/Special/Lab8/Linag/LinearSystem.h"
#include "Labs/Special/Lab8/Base/BoundaryConditions.h"

class PressureCorrectionAssembler {
public:
    static int index(const StaggeredGrid& g, int i, int j) {
        return j * g.pNx() + i;
    }

    static LinearSystem assemble(const StaggeredGrid& g,
                                 const SimulationConfig& cfg,
                                 const UField& uStar,
                                 const VField& vStar,
                                 const ScalarField& dU,
                                 const ScalarField& dV) {
        LinearSystem system;
        system.reset(g.pNx() * g.pNy());

        const int iRef = BoundaryConditions::pressureReferenceI(g);
        const int jRef = BoundaryConditions::pressureReferenceJ(g);

        const double rho = cfg.rho;
        const double dx = g.dx();
        const double dy = g.dy();

        for (int j = 0; j < g.pNy(); ++j) {
            for (int i = 0; i < g.pNx(); ++i) {
                const int row = index(g, i, j);

                if (i == iRef && j == jRef) {
                    system.add(row, row, 1.0);
                    system.add_rhs(row, 0.0);
                    continue;
                }

                const double aW = (i > 0)         ? rho * dU(i, j)     * dy : 0.0;
                const double aE = (i < g.pNx()-1) ? rho * dU(i + 1, j) * dy : 0.0;
                const double aS = (j > 0)         ? rho * dV(i, j)     * dx : 0.0;
                const double aN = (j < g.pNy()-1) ? rho * dV(i, j + 1) * dx : 0.0;
                const double aP = aW + aE + aS + aN;

                const double b = rho * dy * (uStar(i, j) - uStar(i + 1, j))
                               + rho * dx * (vStar(i, j) - vStar(i, j + 1));

                system.add(row, row, aP);
                system.add_rhs(row, b);

                if (i > 0) {
                    system.add(row, index(g, i - 1, j), -aW);
                }
                if (i < g.pNx() - 1) {
                    system.add(row, index(g, i + 1, j), -aE);
                }
                if (j > 0) {
                    system.add(row, index(g, i, j - 1), -aS);
                }
                if (j < g.pNy() - 1) {
                    system.add(row, index(g, i, j + 1), -aN);
                }
            }
        }

        return system;
    }
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_PRESSURECORRECTIONASSEMBLER_H