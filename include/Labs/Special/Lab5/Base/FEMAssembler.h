#ifndef NUMERICAL_METHODS_IN_PHYSICS_FEMASSEMBLER_H
#define NUMERICAL_METHODS_IN_PHYSICS_FEMASSEMBLER_H
#pragma once

#include "Mesh.h"
#include "CylinderHeatProblem.h"
#include "ElementMatrices.h"

#include <Eigen/Sparse>
#include <vector>
#include <cmath>

class FEMAssembler {
public:
    FEMAssembler(const Mesh& mesh, const CylinderHeatProblem& prob)
        : mesh_(mesh), prob_(prob) {}

    void assemble(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const {
        const std::size_t N = mesh_.nodes_.size();

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(mesh_.elements_.size() * 9);

        F = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(N));

        for (const auto& elem : mesh_.elements_) {
            Eigen::Matrix3d Ke;
            Eigen::Vector3d Fe;
            compute_element_matrices(mesh_, prob_, elem, Ke, Fe);

            for (int i = 0; i < 3; ++i) {
                const std::size_t gi = elem.n[i];
                F(static_cast<Eigen::Index>(gi)) += Fe(i);

                for (int j = 0; j < 3; ++j) {
                    const std::size_t gj = elem.n[j];
                    triplets.emplace_back(
                        static_cast<Eigen::Index>(gi),
                        static_cast<Eigen::Index>(gj),
                        Ke(i, j)
                    );
                }
            }
        }

        add_neumann_bottom(F);

        K.resize(static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N));
        K.setFromTriplets(triplets.begin(), triplets.end());
        K.makeCompressed();

        apply_dirichlet(K, F);
    }

private:
    const Mesh& mesh_;
    const CylinderHeatProblem& prob_;

    void add_neumann_bottom(Eigen::VectorXd& F) const {
        const auto& nodes = mesh_.nodes_;
        const auto& ids   = mesh_.bottom_edge_nodes_;

        if (ids.size() < 2 || std::abs(prob_.q) < 1e-14)
            return;

        for (std::size_t k = 0; k + 1 < ids.size(); ++k) {
            const std::size_t i0 = ids[k];
            const std::size_t i1 = ids[k + 1];

            const double r0 = nodes[i0].r;
            const double r1 = nodes[i1].r;
            const double r_edge = 0.5 * (r0 + r1);
            const double length = std::abs(r1 - r0);

            const double contrib = prob_.q * r_edge * length / 2.0;

            F(static_cast<Eigen::Index>(i0)) += contrib;
            F(static_cast<Eigen::Index>(i1)) += contrib;
        }
    }

    void apply_dirichlet(Eigen::SparseMatrix<double>& K,
                         Eigen::VectorXd& F) const
    {
        const Eigen::Index N = K.rows();
        const double H = mesh_.nodes_.back().z;

        std::vector<char> is_dirichlet(static_cast<std::size_t>(N), 0);
        std::vector<double> dirichlet_value(static_cast<std::size_t>(N), 0.0);

        for (Eigen::Index i = 0; i < N; ++i) {
            const Node& node = mesh_.nodes_[static_cast<std::size_t>(i)];

            const bool on_r1  = std::abs(node.r - 1.0) < 1e-14;
            const bool on_top = std::abs(node.z - H)   < 1e-14;

            if (on_r1) {
                is_dirichlet[static_cast<std::size_t>(i)] = 1;
                dirichlet_value[static_cast<std::size_t>(i)] = 0.0;
            } else if (on_top) {
                is_dirichlet[static_cast<std::size_t>(i)] = 1;
                dirichlet_value[static_cast<std::size_t>(i)] = prob_.theta;
            }
        }

        for (int outer = 0; outer < K.outerSize(); ++outer) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(K, outer); it; ++it) {
                const Eigen::Index row = it.row();
                const Eigen::Index col = it.col();

                if (!is_dirichlet[static_cast<std::size_t>(row)] &&
                     is_dirichlet[static_cast<std::size_t>(col)]) {
                    F(row) -= it.value() *
                              dirichlet_value[static_cast<std::size_t>(col)];
                }
            }
        }

        for (int outer = 0; outer < K.outerSize(); ++outer) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(K, outer); it; ++it) {
                const Eigen::Index row = it.row();
                const Eigen::Index col = it.col();

                if (is_dirichlet[static_cast<std::size_t>(row)] ||
                    is_dirichlet[static_cast<std::size_t>(col)]) {
                    it.valueRef() = 0.0;
                }
            }
        }

        for (Eigen::Index i = 0; i < N; ++i) {
            if (!is_dirichlet[static_cast<std::size_t>(i)])
                continue;

            K.coeffRef(i, i) = 1.0;
            F(i) = dirichlet_value[static_cast<std::size_t>(i)];
        }

        K.prune(0.0);
        K.makeCompressed();
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_FEMASSEMBLER_H
