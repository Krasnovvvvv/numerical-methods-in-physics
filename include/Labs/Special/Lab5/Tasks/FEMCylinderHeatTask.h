#ifndef NUMERICAL_METHODS_IN_PHYSICS_FEMCYLINDERHEATTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_FEMCYLINDERHEATTASK_H
#pragma once

#include "Labs/Special/Lab5/Base/CylinderHeatProblem.h"
#include "Labs/Special/Lab5/Base/Mesh.h"
#include "Labs/Special/Lab5/Base/FEMAssembler.h"
#include "Labs/Special/Lab5/Solver/EigenCGSolver.h"

#include <Eigen/Sparse>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <matplot/matplot.h>

class FEMCylinderHeatTask {
public:
    FEMCylinderHeatTask(CylinderHeatProblem problem,
                        EigenCGSolver solver)
        : problem_(std::move(problem)),
          solver_(solver) {}

    void solveSingle(std::size_t nr, std::size_t nz) {
        // Чтобы потом отрисовать nr x nz точек как в МКР,
        // берём на 1 узел больше по каждому направлению
        Mesh mesh;
        mesh.build(problem_.H, nr + 1, nz + 1);

        FEMAssembler assembler(mesh, problem_);

        Eigen::SparseMatrix<double> K;
        Eigen::VectorXd F;
        assembler.assemble(K, F);

        std::cout << "FEM unknowns = " << K.rows() << "\n";

        SolveResult result = solver_.solve(K, F);

        std::cout << "iterations = " << result.iterations << "\n";
        std::cout << "estimated error = "
                  << std::scientific << result.estimatedError << "\n";
        std::cout << "relative residual = "
                  << std::scientific << result.relativeResidual << "\n";

        plotTemperatureField(mesh, result.solution, nr, nz);
    }

private:
    CylinderHeatProblem problem_;
    EigenCGSolver solver_;

    double interpolateAt(const Mesh& mesh,
                         const Eigen::VectorXd& solution,
                         double r,
                         double z) const
    {
        const std::size_t Nr = mesh.nr_;
        const std::size_t Nz = mesh.nz_;

        const double dr = 1.0 / static_cast<double>(Nr - 1);
        const double dz = problem_.H / static_cast<double>(Nz - 1);

        r = std::clamp(r, 0.0, 1.0);
        z = std::clamp(z, 0.0, problem_.H);

        std::size_t i = std::min(
            static_cast<std::size_t>(std::floor(r / dr)),
            Nr - 2
        );
        std::size_t j = std::min(
            static_cast<std::size_t>(std::floor(z / dz)),
            Nz - 2
        );

        const double r0 = mesh.nodes_[mesh.index(i, j)].r;
        const double z0 = mesh.nodes_[mesh.index(i, j)].z;

        double xi  = (r - r0) / dr;
        double eta = (z - z0) / dz;

        xi  = std::clamp(xi,  0.0, 1.0);
        eta = std::clamp(eta, 0.0, 1.0);

        const std::size_t n00 = mesh.index(i,     j);
        const std::size_t n10 = mesh.index(i + 1, j);
        const std::size_t n01 = mesh.index(i,     j + 1);
        const std::size_t n11 = mesh.index(i + 1, j + 1);

        const double T00 = solution(static_cast<Eigen::Index>(n00));
        const double T10 = solution(static_cast<Eigen::Index>(n10));
        const double T01 = solution(static_cast<Eigen::Index>(n01));
        const double T11 = solution(static_cast<Eigen::Index>(n11));

        if (eta <= xi) {
            // Треугольник (n00, n10, n11)
            const double l1 = 1.0 - xi;
            const double l2 = xi - eta;
            const double l3 = eta;
            return l1 * T00 + l2 * T10 + l3 * T11;
        } else {
            // Треугольник (n00, n11, n01)
            const double l1 = 1.0 - eta;
            const double l2 = xi;
            const double l3 = eta - xi;
            return l1 * T00 + l2 * T11 + l3 * T01;
        }
    }

    void plotTemperatureField(const Mesh& mesh,
                              const Eigen::VectorXd& solution,
                              std::size_t nr_mkr,
                              std::size_t nz_mkr)
    {
        using namespace matplot;

        // rCoord(i) = (i + 0.5) * hr,  hr = 1 / (nr - 0.5)
        // zCoord(j) = (j + 0.5) * hz,  hz = H / (nz - 0.5)
        const double hr = 1.0 / (static_cast<double>(nr_mkr) - 0.5);
        const double hz = problem_.H / (static_cast<double>(nz_mkr) - 0.5);

        std::vector<double> r(nr_mkr);
        std::vector<double> z(nz_mkr);

        for (std::size_t i = 0; i < nr_mkr; ++i)
            r[i] = (static_cast<double>(i) + 0.5) * hr;

        for (std::size_t j = 0; j < nz_mkr; ++j)
            z[j] = (static_cast<double>(j) + 0.5) * hz;

        std::vector<std::vector<double>> Tmat(
            nz_mkr, std::vector<double>(nr_mkr, 0.0)
        );

        for (std::size_t j = 0; j < nz_mkr; ++j) {
            for (std::size_t i = 0; i < nr_mkr; ++i) {
                Tmat[j][i] = interpolateAt(mesh, solution, r[i], z[j]);
            }
        }

        auto [R, Z] = meshgrid(r, z);

        figure();
        auto ax = gca();
        contourf(R, Z, Tmat, 40);
        colorbar();

        ax->title("FEM Temperature field T(r,z)");
        ax->xlabel("r");
        ax->ylabel("z");
        ax->xlim({r.front(), r.back()});
        ax->ylim({z.front(), z.back()});
        ax->grid(false);

        show();
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_FEMCYLINDERHEATTASK_H
