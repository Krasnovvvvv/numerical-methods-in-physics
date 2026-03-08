#ifndef NUMERICAL_METHODS_IN_PHYSICS_CYLINDERHEATGENERATOR_H
#define NUMERICAL_METHODS_IN_PHYSICS_CYLINDERHEATGENERATOR_H

#pragma once
#include "CylinderHeatProblem.h"

#include <Eigen/Sparse>
#include <vector>
#include <utility>
#include <stdexcept>

class CylinderHeatGenerator {
public:
    explicit CylinderHeatGenerator(CylinderHeatProblem problem)
        : problem_(std::move(problem))
    {
        if (problem_.nr == 0 || problem_.nz == 0)
            throw std::runtime_error("nr and nz must be positive.");
        if (problem_.H <= 0.0)
            throw std::runtime_error("H must be positive.");
    }

    [[nodiscard]] std::pair<
        Eigen::SparseMatrix<double>,
        Eigen::VectorXd> generateSLAE() const {
        const size_t N = unknowns();
        const double dr = hr();
        const double dz = hz();
        const double inv_dr2 = 1.0 / (dr * dr);
        const double inv_dz2 = 1.0 / (dz * dz);

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(5 * N);

        Eigen::VectorXd b = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(N));

        for (size_t j = 0; j < problem_.nz; ++j) {
            const double z = zCoord(j);

            for (size_t i = 0; i < problem_.nr; ++i) {
                const double r  = rCoord(i);
                const double rw = static_cast<double>(i) * dr;         // r_{i-1/2}
                const double re = (static_cast<double>(i) + 1.0) * dr; // r_{i+1/2}

                const size_t k = index(i, j);
                double diag = 0.0;
                double rhs = r * problem_.rhs(r, z);

                if (i > 0) {
                    const double aw = rw * inv_dr2;
                    diag += aw;
                    triplets.emplace_back(static_cast<Eigen::Index>(k),
                                          static_cast<Eigen::Index>(index(i - 1, j)),
                                          -aw);
                }

                if (i + 1 < problem_.nr) {
                    const double ae = re * inv_dr2;
                    diag += ae;
                    triplets.emplace_back(static_cast<Eigen::Index>(k),
                                          static_cast<Eigen::Index>(index(i + 1, j)),
                                          -ae);
                } else {
                    // r = 1, T = 0
                    // T_g = -T_i
                    const double ae = re * inv_dr2;
                    diag += 2.0 * ae;
                }

                if (j > 0) {
                    const double as = r * inv_dz2;
                    diag += as;
                    triplets.emplace_back(static_cast<Eigen::Index>(k),
                                          static_cast<Eigen::Index>(index(i, j - 1)),
                                          -as);
                } else {
                    // z = 0, dT/dz = -q
                    rhs += r * problem_.q / dz;
                }

                if (j + 1 < problem_.nz) {
                    const double an = r * inv_dz2;
                    diag += an;
                    triplets.emplace_back(static_cast<Eigen::Index>(k),
                                          static_cast<Eigen::Index>(index(i, j + 1)),
                                          -an);
                } else {
                    // z = H, T = theta
                    // T_g = 2*theta - T_i
                    diag += 2.0 * r * inv_dz2;
                    rhs += 2.0 * r * problem_.theta * inv_dz2;
                }

                triplets.emplace_back(static_cast<Eigen::Index>(k),
                                      static_cast<Eigen::Index>(k),
                                      diag);

                b(static_cast<Eigen::Index>(k)) = rhs;
            }
        }

        Eigen::SparseMatrix<double> A(static_cast<Eigen::Index>(N),
                                      static_cast<Eigen::Index>(N));
        A.setFromTriplets(triplets.begin(), triplets.end());
        A.makeCompressed();

        return {A, b};
    }

    [[nodiscard]] std::vector<double> radialPoints() const {
        std::vector<double> r(problem_.nr);
        for (size_t i = 0; i < problem_.nr; ++i)
            r[i] = rCoord(i);
        return r;
    }

    [[nodiscard]] std::vector<double> axialPoints() const {
        std::vector<double> z(problem_.nz);
        for (size_t j = 0; j < problem_.nz; ++j)
            z[j] = zCoord(j);
        return z;
    }

    [[nodiscard]] size_t nr() const { return problem_.nr; }
    [[nodiscard]] size_t nz() const { return problem_.nz; }
    [[nodiscard]] double H() const { return problem_.H; }

    [[nodiscard]] size_t index(size_t i, size_t j) const {
        return j * problem_.nr + i;
    }

    [[nodiscard]] double rCoord(size_t i) const {
        return (static_cast<double>(i) + 0.5) * hr();
    }

    [[nodiscard]] double zCoord(size_t j) const {
        return (static_cast<double>(j) + 0.5) * hz();
    }

private:
    CylinderHeatProblem problem_;

    [[nodiscard]] size_t unknowns() const {
        return problem_.nr * problem_.nz;
    }

    [[nodiscard]] double hr() const {
        return 1.0 / (static_cast<double>(problem_.nr) - 0.5);
    }

    [[nodiscard]] double hz() const {
        return problem_.H / (static_cast<double>(problem_.nz) - 0.5);
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_CYLINDERHEATGENERATOR_H