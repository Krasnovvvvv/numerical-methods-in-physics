#ifndef NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONDIFFUSIONSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONDIFFUSIONSOLVER_H
#pragma once

#include <stdexcept>

#include "Labs/Special/Lab7/Base/ConvectionDiffusionAssembler.h"
#include "IterativeSolver.h"
#include "Labs/Special/Lab7/Base/ScalarField.h"

struct NonlinearSolveResult {
    int outerIterations = 0;
    double maxDelta = 0.0;
    bool converged = false;
    SolverResult linearResult;
    std::vector<double> history;
};

class ConvectionDiffusionSolver {
public:
    explicit ConvectionDiffusionSolver(const IterativeSolver& linearSolver_,
                                       int maxOuterIterations_ = 100,
                                       double outerTolerance_ = 1e-10)
        : linearSolver(linearSolver_),
          maxOuterIterations(maxOuterIterations_),
          outerTolerance(outerTolerance_) {}

    NonlinearSolveResult solve(const ConvectionDiffusionProblem& pr, ScalarField& phi) const {
        validateCompatible(pr.mesh, phi.mesh);

        NonlinearSolveResult result;

        if (pr.scheme == ConvectionScheme::Upwind) {
            ScalarField oldPhi = phi;
            LinearSystem sys = ConvectionDiffusionAssembler::build(pr, oldPhi);
            result.linearResult = linearSolver.solve(sys, phi);
            result.maxDelta = phi.maxAbsDiff(oldPhi);
            result.history.push_back(result.maxDelta);
            result.outerIterations = 1;
            result.converged = result.linearResult.converged;
            return result;
        }

        for (int iter = 1; iter <= maxOuterIterations; ++iter) {
            ScalarField phiLag = phi;
            LinearSystem sys = ConvectionDiffusionAssembler::build(pr, phiLag);
            result.linearResult = linearSolver.solve(sys, phi);
            result.maxDelta = phi.maxAbsDiff(phiLag);
            result.history.push_back(result.maxDelta);
            result.outerIterations = iter;
            result.converged = result.linearResult.converged && (result.maxDelta <= outerTolerance);
            if (result.converged) return result;
        }

        return result;
    }

private:
    IterativeSolver linearSolver;
    int maxOuterIterations;
    double outerTolerance;

    static void validateCompatible(const Mesh& a, const Mesh& b) {
        if (a.Nx != b.Nx || a.Ny != b.Ny || a.Lx != b.Lx || a.Ly != b.Ly) {
            throw std::runtime_error("ConvectionDiffusionSolver: incompatible meshes.");
        }
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONDIFFUSIONSOLVER_H