#ifndef NUMERICAL_METHODS_IN_PHYSICS_BOUNDARYCONDITION_H
#define NUMERICAL_METHODS_IN_PHYSICS_BOUNDARYCONDITION_H
#pragma once
#include <algorithm>
#include "Labs/Special/Lab8/Base/StaggeredGrid.h"
#include "Labs/Special/Lab8/Base/Config/SimConfig.h"
#include "Labs/Special/Lab8/Base/Field/ScalarField.h"
#include "Labs/Special/Lab8/Base/Field/UField.h"
#include "Labs/Special/Lab8/Base/Field/VField.h"

struct BoundaryConditions {
    static double inletU(double y, const SimulationConfig& cfg) {
        if (cfg.inletType == InletType::Uniform) {
            return cfg.uMean;
        }
        const double eta = std::clamp(y / cfg.Ly, 0.0, 1.0);
        return 6.0 * cfg.uMean * eta * (1.0 - eta);
    }

    static void applyVelocityBC(const StaggeredGrid& g,
                                const SimulationConfig& cfg,
                                UField& u,
                                VField& v) {
        for (int j = 0; j < g.uNy(); ++j) {
            const double y = (j + 0.5) * g.dy();
            u(0, j) = inletU(y, cfg);
            u(g.uNx() - 1, j) = u(g.uNx() - 2, j);
        }

        for (int i = 0; i < g.vNx(); ++i) {
            v(i, 0) = 0.0;
            v(i, g.vNy() - 1) = 0.0;
        }

        for (int j = 1; j < g.vNy() - 1; ++j) {
            v(0, j) = cfg.vInlet;
            v(g.vNx() - 1, j) = v(g.vNx() - 2, j);
        }
    }

    static int pressureReferenceI(const StaggeredGrid& g) {
        return g.pNx() - 1;
    }

    static int pressureReferenceJ(const StaggeredGrid& g) {
        return g.pNy() / 2;
    }

    static void shiftPressureToOutletReference(const StaggeredGrid& g,
                                               const SimulationConfig& cfg,
                                               ScalarField& p) {
        const int iRef = pressureReferenceI(g);
        const int jRef = pressureReferenceJ(g);
        const double shift = cfg.pOutlet - p(iRef, jRef);
        for (int j = 0; j < g.pNy(); ++j) {
            for (int i = 0; i < g.pNx(); ++i) {
                p(i, j) += shift;
            }
        }
    }
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_BOUNDARYCONDITION_H