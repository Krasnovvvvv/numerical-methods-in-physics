#ifndef NUMERICAL_METHODS_IN_PHYSICS_CYLINDERHEATTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_CYLINDERHEATTASK_H

#pragma once
#include "Labs/Special/Lab5/Base/CylinderHeatProblem.h"
#include "Labs/Special/Lab5/Base/CylinderHeatGenerator.h"
#include "Labs/Special/Lab5/Solver/EigenCGSolver.h"
#include "Helpers/Plotter.h"

#include <Eigen/Sparse>
#include <matplot/matplot.h>
#include <iostream>
#include <string>
#include <vector>


class CylinderHeatTask {
public:
    CylinderHeatTask(CylinderHeatProblem baseProblem,
                     EigenCGSolver solver,
                     Plotter& plotter)
        : baseProblem_(std::move(baseProblem)),
          solver_(solver),
          plotter_(plotter) {}

    void solveSingle(size_t nr, size_t nz)
    {
        CylinderHeatProblem problem = baseProblem_;
        problem.nr = nr;
        problem.nz = nz;

        CylinderHeatGenerator generator(problem);
        auto [A, b] = generator.generateSLAE();
        SolveResult result = solver_.solve(A, b);

        std::cout << "\nSingle run:\n";
        std::cout << "unknowns          = "
        << problem.nr * problem.nz << "\n";
        std::cout << "iterations        = "
        << result.iterations << "\n";
        std::cout << "estimated error   = "
        << std::scientific << result.estimatedError << "\n";
        std::cout << "relative residual = "
        << std::scientific << result.relativeResidual << "\n";

        plotAxisProfile(generator, result.solution);
        plotRadialSlices(generator, result.solution);

        plotTemperatureField(generator, result.solution);

    }

private:
    CylinderHeatProblem baseProblem_;
    EigenCGSolver solver_;
    Plotter& plotter_;

    void plotAxisProfile(const CylinderHeatGenerator& generator,
                         const Eigen::VectorXd& solution)
    {
        std::vector<double> z = generator.axialPoints();
        std::vector<double> T(generator.nz());

        for (size_t j = 0; j < generator.nz(); ++j) {
            T[j] = solution(static_cast<Eigen::Index>
                (generator.index(0, j)));
        }

        plotter_.plot(
            z, T, "T(r≈0, z)",
            "z", "T", false);
    }

    void plotRadialSlices(const CylinderHeatGenerator& generator,
                          const Eigen::VectorXd& solution)
    {
        std::vector<double> r = generator.radialPoints();

        const size_t jBottom = 0;
        const size_t jMiddle = generator.nz() / 2;
        const size_t jTop = generator.nz() - 1;

        std::vector<double> Tbottom(generator.nr());
        std::vector<double> Tmiddle(generator.nr());
        std::vector<double> Ttop(generator.nr());

        for (size_t i = 0; i < generator.nr(); ++i) {
            Tbottom[i] = solution(static_cast<Eigen::Index>
                (generator.index(i, jBottom)));
            Tmiddle[i] = solution(static_cast<Eigen::Index>
                (generator.index(i, jMiddle)));
            Ttop[i]    = solution(static_cast<Eigen::Index>
                (generator.index(i, jTop)));
        }

        std::vector<std::vector<double>> xs = {r, r, r};
        std::vector<std::vector<double>> ys = {Tbottom, Tmiddle, Ttop};
        std::vector<std::string> labels = {
            "z = z_{min}",
            "z = z_{mid}",
            "z = z_{max}"
        };

        plotter_.plot(xs, ys, labels, "r", "T");
    }

    void plotTemperatureField(const CylinderHeatGenerator& generator,
                          const Eigen::VectorXd& solution)
    {
        using namespace matplot;

        size_t nr = generator.nr();
        size_t nz = generator.nz();

        std::vector<double> r = generator.radialPoints(); // r in (0,1)
        std::vector<double> z = generator.axialPoints();  // z in (0,H)

        std::vector Tmat(nz, std::vector<double>(nr));
        for (size_t j = 0; j < nz; ++j) {
            for (size_t i = 0; i < nr; ++i) {
                size_t k = generator.index(i, j);
                Tmat[j][i] = solution(static_cast<Eigen::Index>(k));
            }
        }

        auto [R, Z] = meshgrid(r, z);

        figure();
        auto ax = gca();
        contourf(R, Z, Tmat, 40);
        colorbar();

        ax->title("Temperature field T(r,z)");
        ax->xlabel("r");
        ax->ylabel("z");
        ax->xlim({r.front(), r.back()});
        ax->ylim({z.front(), z.back()});
        ax->grid(true);

        show();
    }

};

#endif // NUMERICAL_METHODS_IN_PHYSICS_CYLINDERHEATTASK_H