#ifndef NUMERICAL_METHODS_IN_PHYSICS_EIGENSPARSESOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_EIGENSPARSESOLVER_H
#pragma once
#include "LinearSolver.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <stdexcept>

class EigenSparseSolver final : public LinearSolver {
public:
    std::vector<double> solve(const LinearSystem& system) override {
        using SparseMatrix = Eigen::SparseMatrix<double>;
        using Triplet = Eigen::Triplet<double>;

        std::vector<Triplet> triplets;
        triplets.reserve(system.vals.size());
        for (size_t k = 0; k < system.vals.size(); ++k) {
            triplets.emplace_back(system.rows[k], system.cols[k], system.vals[k]);
        }

        SparseMatrix A(system.n, system.n);
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::VectorXd b(system.n);
        for (int i = 0; i < system.n; ++i) {
            b[i] = system.rhs[static_cast<size_t>(i)];
        }

        Eigen::SparseLU<SparseMatrix> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("EigenSparseSolver: factorization failed");
        }

        Eigen::VectorXd x = solver.solve(b);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("EigenSparseSolver: solve failed");
        }

        return std::vector<double>(x.data(), x.data() + x.size());
    }
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_EIGENSPARSESOLVER_H