#ifndef NUMERICAL_METHODS_IN_PHYSICS_ITERATIVESOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_ITERATIVESOLVER_H
#pragma once

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>

#include <vector>
#include <stdexcept>
#include <string>
#include <cmath>
#include <algorithm>

#include "Labs/Special/Lab7/Base/LinearSystem.h"
#include "Labs/Special/Lab7/Base/ScalarField.h"

struct SolverResult {
    int iterations = 0;
    double residual = 0.0;
    double maxDelta = 0.0;
    bool converged = false;
    std::string solverName;
};

enum class EigenSolverType {
    SparseLU,
    BiCGSTAB
};

class IterativeSolver {
public:
    explicit IterativeSolver(EigenSolverType type_ = EigenSolverType::SparseLU,
                             int maxIterations_ = 5000,
                             double tolerance_ = 1e-10)
        : type(type_), maxIterations(maxIterations_), tolerance(tolerance_) {}

    SolverResult solve(const LinearSystem& A, ScalarField& phi) const {
        validateCompatible(A.mesh, phi.mesh);

        Eigen::SparseMatrix<double> M = buildMatrix(A);
        Eigen::VectorXd b = buildRhs(A);
        Eigen::VectorXd x = buildInitialGuess(phi);
        Eigen::VectorXd xOld = x;

        SolverResult result;

        if (type == EigenSolverType::SparseLU) {
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.analyzePattern(M);
            solver.factorize(M);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("SparseLU factorization failed.");
            }
            x = solver.solve(b);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("SparseLU solve failed.");
            }
            result.iterations = 1;
            result.converged = true;
            result.solverName = "Eigen::SparseLU";
        } else {
            Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
            solver.setMaxIterations(maxIterations);
            solver.setTolerance(tolerance);
            solver.compute(M);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("BiCGSTAB compute failed.");
            }
            x = solver.solveWithGuess(b, x);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("BiCGSTAB solve failed.");
            }
            result.iterations = solver.iterations();
            result.converged = (solver.error() <= tolerance);
            result.solverName = "Eigen::BiCGSTAB";
        }

        const double bNorm = std::max(b.norm(), 1e-14);
        result.residual = (M * x - b).norm() / bNorm;
        result.maxDelta = (x - xOld).cwiseAbs().maxCoeff();

        writeBack(phi, x);
        return result;
    }

private:
    EigenSolverType type;
    int maxIterations;
    double tolerance;

    static void validateCompatible(const Mesh& a, const Mesh& b) {
        if (a.Nx != b.Nx || a.Ny != b.Ny || a.Lx != b.Lx || a.Ly != b.Ly) {
            throw std::runtime_error("IterativeSolver: incompatible meshes.");
        }
    }

    static Eigen::SparseMatrix<double> buildMatrix(const LinearSystem& A) {
        const Mesh& m = A.mesh;
        const int n = m.cellCount();
        const double eps = 1e-14;

        Eigen::SparseMatrix<double> M(n, n);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(static_cast<size_t>(5 * n));

        for (int j = 0; j < m.Ny; ++j) {
            for (int i = 0; i < m.Nx; ++i) {
                const int p = m.id(i, j);
                const size_t ps = static_cast<size_t>(p);

                if (std::abs(A.aP[ps]) < eps) {
                    throw std::runtime_error("Zero diagonal coefficient encountered.");
                }

                triplets.emplace_back(p, p, A.aP[ps]);

                if (m.hasWest(i)  && std::abs(A.aW[ps]) > eps) triplets.emplace_back(p, m.id(i - 1, j), -A.aW[ps]);
                if (m.hasEast(i)  && std::abs(A.aE[ps]) > eps) triplets.emplace_back(p, m.id(i + 1, j), -A.aE[ps]);
                if (m.hasSouth(j) && std::abs(A.aS[ps]) > eps) triplets.emplace_back(p, m.id(i, j - 1), -A.aS[ps]);
                if (m.hasNorth(j) && std::abs(A.aN[ps]) > eps) triplets.emplace_back(p, m.id(i, j + 1), -A.aN[ps]);
            }
        }

        M.setFromTriplets(triplets.begin(), triplets.end());
        M.makeCompressed();
        return M;
    }

    static Eigen::VectorXd buildRhs(const LinearSystem& A) {
        Eigen::VectorXd b(A.mesh.cellCount());
        for (int k = 0; k < A.mesh.cellCount(); ++k) b[k] = A.b[static_cast<size_t>(k)];
        return b;
    }

    static Eigen::VectorXd buildInitialGuess(const ScalarField& phi) {
        Eigen::VectorXd x(phi.size());
        for (int k = 0; k < phi.size(); ++k) x[k] = phi.data[static_cast<size_t>(k)];
        return x;
    }

    static void writeBack(ScalarField& phi, const Eigen::VectorXd& x) {
        for (int k = 0; k < phi.size(); ++k) phi.data[static_cast<size_t>(k)] = x[k];
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_ITERATIVESOLVER_H