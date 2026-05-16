#ifndef NUMERICAL_METHODS_IN_PHYSICS_ADVECTIONDIFFUSIONTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_ADVECTIONDIFFUSIONTASK_H
#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <optional>

#include "Labs/Special/Lab7/Base/Mesh.h"
#include "Labs/Special/Lab7/Base/ScalarField.h"
#include "Labs/Special/Lab7/Base/BoundaryCondition.h"
#include "Labs/Special/Lab7/Base/ConvectionDiffusionProblem.h"
#include "Labs/Special/Lab7/Solver/ConvectionDiffusionSolver.h"
#include "Labs/Special/Lab7/Solver/IterativeSolver.h"
#include "Labs/Special/Lab7/Base/CFDPlotter.h"

struct AdvectionDiffusionTaskConfig {
    int Nx = 120;
    int Ny = 120;
    double Lx = 1.0;
    double Ly = 1.0;

    double rho = 1.0;
    double u = 2.0;
    double v = 2.0;
    double Gamma = 0.0;

    double Su = 0.0;
    double Sp = 0.0;

    BoundaryCondition left   = BoundaryCondition::dirichlet(1.0);
    BoundaryCondition right  = BoundaryCondition::outlet();
    BoundaryCondition bottom = BoundaryCondition::dirichlet(0.0);
    BoundaryCondition top    = BoundaryCondition::outlet();

    double initialValue = 0.0;

    ConvectionScheme scheme = ConvectionScheme::TVD;
    TVDLimiterType limiter  = TVDLimiterType::VanLeer;
    double swebyBeta        = 1.5;

    EigenSolverType linearSolverType = EigenSolverType::SparseLU;
    int linearMaxIterations = 5000;
    double linearTolerance = 1e-12;

    int outerMaxIterations = 100;
    double outerTolerance = 1e-10;

    bool showHeatmap = true;
    bool showProfileX = false;
    bool showProfileY = false;
    bool showHistory = false;

    std::optional<double> profileY;
    std::optional<double> profileX;

    std::string title = "Advection-Diffusion Task";
};

struct AdvectionDiffusionTaskResult {
    ScalarField field;
    NonlinearSolveResult solve;
};

class AdvectionDiffusionTask {
public:
    static AdvectionDiffusionTaskResult solve(const AdvectionDiffusionTaskConfig& cfg) {
        Mesh mesh(cfg.Nx, cfg.Ny, cfg.Lx, cfg.Ly);

        ConvectionDiffusionProblem pr(mesh);
        pr.rho = cfg.rho;
        pr.u = cfg.u;
        pr.v = cfg.v;
        pr.Gamma = cfg.Gamma;
        pr.Su = cfg.Su;
        pr.Sp = cfg.Sp;
        pr.left = cfg.left;
        pr.right = cfg.right;
        pr.bottom = cfg.bottom;
        pr.top = cfg.top;
        pr.scheme = cfg.scheme;
        pr.limiter = cfg.limiter;
        pr.swebyBeta = cfg.swebyBeta;

        IterativeSolver linearSolver(
            cfg.linearSolverType,
            cfg.linearMaxIterations,
            cfg.linearTolerance
        );

        ConvectionDiffusionSolver solver(
            linearSolver,
            cfg.outerMaxIterations,
            cfg.outerTolerance
        );

        ScalarField phi(mesh, cfg.initialValue);
        NonlinearSolveResult nonlinear = solver.solve(pr, phi);

        std::cout << "\n=== " << cfg.title << " ===\n";
        std::cout << "Scheme          : " << schemeName(cfg.scheme, cfg.limiter) << '\n';
        std::cout << "Gamma           : " << cfg.Gamma << '\n';
        std::cout << "Outer iterations: " << nonlinear.outerIterations << '\n';
        std::cout << "Outer max delta : " << nonlinear.maxDelta << '\n';
        std::cout << "Linear residual : " << nonlinear.linearResult.residual << '\n';
        std::cout << "Linear solver   : " << nonlinear.linearResult.solverName << '\n';
        std::cout << "Converged       : " << (nonlinear.converged ? "yes" : "no") << '\n';

        return {phi, nonlinear};
    }

    static void plotStandard(const AdvectionDiffusionTaskConfig& cfg,
                             const AdvectionDiffusionTaskResult& result) {
        const double py = cfg.profileY.value_or(0.5 * cfg.Ly);
        const double px = cfg.profileX.value_or(0.5 * cfg.Lx);

        if (cfg.showHeatmap) {
            CFDPlotter::plotHeatmap(result.field, cfg.title);
        }
        if (cfg.showProfileX) {
            CFDPlotter::plotProfileX(
                result.field,
                py,
                "main",
                cfg.title + " : c(x, y = const)"
            );
        }
        if (cfg.showProfileY) {
            CFDPlotter::plotProfileY(
                result.field,
                px,
                "main",
                cfg.title + " : c(x = const, y)"
            );
        }
        if (cfg.showHistory) {
            CFDPlotter::plotHistory(
                result.solve.history,
                cfg.title + " : outer iteration history"
            );
        }
    }

    static AdvectionDiffusionTaskResult run(const AdvectionDiffusionTaskConfig& cfg) {
        auto result = solve(cfg);
        plotStandard(cfg, result);
        return result;
    }

    static std::string limiterName(TVDLimiterType lim) {
        switch (lim) {
            case TVDLimiterType::VanLeer:   return "VanLeer";
            case TVDLimiterType::VanAlbada: return "VanAlbada";
            case TVDLimiterType::MinMod:    return "MinMod";
            case TVDLimiterType::Superbee:  return "Superbee";
            case TVDLimiterType::Sweby:     return "Sweby";
            case TVDLimiterType::QUICK:     return "QUICK";
            case TVDLimiterType::UMIST:     return "UMIST";
            default:                        return "Unknown";
        }
    }

    static std::string schemeName(ConvectionScheme scheme, TVDLimiterType lim) {
        if (scheme == ConvectionScheme::Upwind) {
            return "Upwind";
        }
        return "TVD-" + limiterName(lim);
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_ADVECTIONDIFFUSIONTASK_H