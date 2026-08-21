#ifndef NUMERICAL_METHODS_IN_PHYSICS_SIMPLECSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_SIMPLECSOLVER_H
#pragma once
#include <cmath>
#include <iostream>
#include <algorithm>
#include "Labs/Special/Lab8/Base/StaggeredGrid.h"
#include "Labs/Special/Lab8/Base/Config/SimConfig.h"
#include "Labs/Special/Lab8/Base/Field/ScalarField.h"
#include "Labs/Special/Lab8/Base/Field/UField.h"
#include "Labs/Special/Lab8/Base/Field/VField.h"
#include "Labs/Special/Lab8/Linag/LinearSolver.h"
#include "Labs/Special/Lab8/Base/BoundaryConditions.h"
#include "Labs/Special/Lab8/FVM/MomentumUAssembler.h"
#include "Labs/Special/Lab8/FVM/MomentumVAssembler.h"
#include "Labs/Special/Lab8/FVM/PressureCorrectionAssembler.h"
#include "FlowResult.h"

class SimpleCSolver {
public:
    SimpleCSolver(StaggeredGrid grid, SimulationConfig cfg, LinearSolver& linearSolver)
        : grid_(grid), cfg_(cfg), linearSolver_(linearSolver),
          p_(grid.pNx(), grid.pNy(), 0.0),
          pCorr_(grid.pNx(), grid.pNy(), 0.0),
          u_(grid.uNx(), grid.uNy(), cfg.uMean),
          v_(grid.vNx(), grid.vNy(), 0.0),
          uStar_(grid.uNx(), grid.uNy(), cfg.uMean),
          vStar_(grid.vNx(), grid.vNy(), 0.0),
          dU_(grid.uNx(), grid.uNy(), 0.0),
          dV_(grid.vNx(), grid.vNy(), 0.0) {}

    FlowResult solve() {
        initialize();
        Residuals current;
        int iter = 0;

        for (iter = 1; iter <= cfg_.maxIterations; ++iter) {
            const UField uOld = u_;
            const VField vOld = v_;
            const ScalarField pOld = p_;

            BoundaryConditions::applyVelocityBC(grid_, cfg_, u_, v_);

            solveMomentumU();
            solveMomentumV();
            solvePressureCorrection();
            correctPressureVelocity();

            current = computeResiduals(uOld, vOld, pOld);
            if (iter % cfg_.printEvery == 0 || iter == 1) {
                std::cout << "iter = " << iter
                          << ", mass = " << current.mass
                          << ", du = " << current.du
                          << ", dv = " << current.dv
                          << ", dp = " << current.dp
                          << std::endl;
            }

            if (current.mass < cfg_.tolMass &&
                current.du < cfg_.tolU &&
                current.dv < cfg_.tolV &&
                current.dp < cfg_.tolP) {
                break;
            }
        }

        FlowResult result{grid_, cfg_, p_, u_, v_, current, iter, computeMaxYPlus()};
        return result;
    }

private:
    void initialize() {
        p_.fill(cfg_.pOutlet);
        u_.fill(cfg_.uMean);
        v_.fill(0.0);
        uStar_.fill(cfg_.uMean);
        vStar_.fill(0.0);

        BoundaryConditions::applyVelocityBC(grid_, cfg_, u_, v_);
        BoundaryConditions::shiftPressureToOutletReference(grid_, cfg_, p_);
    }

    void solveMomentumU() {
        auto assembled = MomentumUAssembler::assemble(grid_, cfg_, u_, v_, p_);
        dU_ = assembled.d;
        const auto x = linearSolver_.solve(assembled.system);

        uStar_ = u_;
        for (int j = 0; j < grid_.uNy(); ++j) {
            for (int i = 1; i < grid_.uNx() - 1; ++i) {
                uStar_(i, j) = x[static_cast<size_t>(MomentumUAssembler::index(grid_, i, j))];
            }
        }
        BoundaryConditions::applyVelocityBC(grid_, cfg_, uStar_, vStar_);
    }

    void solveMomentumV() {
        auto assembled = MomentumVAssembler::assemble(grid_, cfg_, uStar_, v_, p_);
        dV_ = assembled.d;
        const auto x = linearSolver_.solve(assembled.system);

        vStar_ = v_;
        for (int j = 1; j < grid_.vNy() - 1; ++j) {
            for (int i = 0; i < grid_.vNx(); ++i) {
                vStar_(i, j) = x[static_cast<size_t>(MomentumVAssembler::index(grid_, i, j))];
            }
        }
        BoundaryConditions::applyVelocityBC(grid_, cfg_, uStar_, vStar_);
    }

    void solvePressureCorrection() {
        auto system = PressureCorrectionAssembler::assemble(grid_, cfg_, uStar_, vStar_, dU_, dV_);
        const auto x = linearSolver_.solve(system);
        for (int j = 0; j < grid_.pNy(); ++j) {
            for (int i = 0; i < grid_.pNx(); ++i) {
                pCorr_(i, j) = x[static_cast<size_t>(PressureCorrectionAssembler::index(grid_, i, j))];
            }
        }
    }

    void correctPressureVelocity() {
        for (int j = 0; j < grid_.pNy(); ++j) {
            for (int i = 0; i < grid_.pNx(); ++i) {
                p_(i, j) += cfg_.alphaP * pCorr_(i, j);
            }
        }
        BoundaryConditions::shiftPressureToOutletReference(grid_, cfg_, p_);

        u_ = uStar_;
        for (int j = 0; j < grid_.uNy(); ++j) {
            for (int i = 1; i < grid_.uNx() - 1; ++i) {
                u_(i, j) = uStar_(i, j) + dU_(i, j) * (pCorr_(i - 1, j) - pCorr_(i, j));
            }
        }

        v_ = vStar_;
        for (int j = 1; j < grid_.vNy() - 1; ++j) {
            for (int i = 0; i < grid_.vNx(); ++i) {
                v_(i, j) = vStar_(i, j) + dV_(i, j) * (pCorr_(i, j - 1) - pCorr_(i, j));
            }
        }

        BoundaryConditions::applyVelocityBC(grid_, cfg_, u_, v_);
    }

    Residuals computeResiduals(const UField& uOld,
                               const VField& vOld,
                               const ScalarField& pOld) const {
        Residuals res;
        for (int j = 0; j < grid_.pNy(); ++j) {
            for (int i = 0; i < grid_.pNx(); ++i) {
                const double m = cfg_.rho * grid_.dy() * (u_(i, j) - u_(i + 1, j))
                               + cfg_.rho * grid_.dx() * (v_(i, j) - v_(i, j + 1));
                res.mass = std::max(res.mass, std::abs(m));
            }
        }
        for (int j = 0; j < grid_.uNy(); ++j) {
            for (int i = 0; i < grid_.uNx(); ++i) {
                res.du = std::max(res.du, std::abs(u_(i, j) - uOld(i, j)));
            }
        }
        for (int j = 0; j < grid_.vNy(); ++j) {
            for (int i = 0; i < grid_.vNx(); ++i) {
                res.dv = std::max(res.dv, std::abs(v_(i, j) - vOld(i, j)));
            }
        }
        for (int j = 0; j < grid_.pNy(); ++j) {
            for (int i = 0; i < grid_.pNx(); ++i) {
                res.dp = std::max(res.dp, std::abs(p_(i, j) - pOld(i, j)));
            }
        }
        return res;
    }

    double computeMaxYPlus() const {
        const double mu = cfg_.mu();
        const double yP = 0.5 * grid_.dy();
        double maxYPlus = 0.0;

        for (int i = 0; i < grid_.pNx(); ++i) {
            const double uBottom = 0.5 * (u_(i, 0) + u_(i + 1, 0));
            const double tauBottom = mu * std::abs(uBottom) / yP;
            const double uTauBottom = std::sqrt(tauBottom / cfg_.rho);
            maxYPlus = std::max(maxYPlus, cfg_.rho * uTauBottom * yP / mu);

            const int jTop = grid_.pNy() - 1;
            const double uTop = 0.5 * (u_(i, jTop) + u_(i + 1, jTop));
            const double tauTop = mu * std::abs(uTop) / yP;
            const double uTauTop = std::sqrt(tauTop / cfg_.rho);
            maxYPlus = std::max(maxYPlus, cfg_.rho * uTauTop * yP / mu);
        }
        return maxYPlus;
    }

private:
    StaggeredGrid grid_;
    SimulationConfig cfg_;
    LinearSolver& linearSolver_;
    ScalarField p_;
    ScalarField pCorr_;
    UField u_;
    VField v_;
    UField uStar_;
    VField vStar_;
    ScalarField dU_;
    ScalarField dV_;
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_SIMPLECSOLVER_H