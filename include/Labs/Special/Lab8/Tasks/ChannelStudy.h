#ifndef NUMERICAL_METHODS_IN_PHYSICS_CHANNELSTUDY_H
#define NUMERICAL_METHODS_IN_PHYSICS_CHANNELSTUDY_H
#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include "Labs/Special/Lab8/Base/Config/CaseSetup.h"
#include "Labs/Special/Lab8/Base/StaggeredGrid.h"
#include "Labs/Special/Lab8/Linag/EigenSparseSolver.h"
#include "Labs/Special/Lab8/Solver/SIMPLECSolver.h"
#include "Labs/Special/Lab8/Base/Post/Plotter.h"

class ChannelStudy {
public:
    static FlowResult solveCase(const SimulationConfig& cfg) {
        StaggeredGrid grid(cfg.Nx, cfg.Ny, cfg.Lx, cfg.Ly);
        EigenSparseSolver linearSolver;
        SimpleCSolver solver(grid, cfg, linearSolver);
        return solver.solve();
    }

    static void runAll() {
        std::cout << std::fixed << std::setprecision(6);

        auto base = make_base_case();
        auto uniformCfg = base;
        uniformCfg.inletType = InletType::Uniform;

        auto profileCfg = base;
        profileCfg.inletType = InletType::Parabolic;

        /*std::cout << "=== Base low-Re run: uniform inlet ===" << std::endl;
        auto uniformResult = solveCase(uniformCfg);
        std::cout << "iterations = " << uniformResult.iterations
                  << ", max y+ = " << uniformResult.maxYPlus << std::endl;

        std::cout << "=== Base low-Re run: parabolic inlet ===" << std::endl;
        auto profileResult = solveCase(profileCfg);
        std::cout << "iterations = " << profileResult.iterations
                  << ", max y+ = " << profileResult.maxYPlus << std::endl;

        //Plotter::plotFieldOverview(uniformResult, "Low Re, uniform inlet");
        Plotter::plotFieldOverview(profileResult, "Low Re, profile inlet");
        Plotter::plotInletComparison(uniformResult, profileResult);

        std::vector<double> lengths = {0.5, 1.0, 4.0};
        std::vector<FlowResult> outletResults;
        std::cout << "=== Outlet study ===" << std::endl;
        for (double L : lengths) {
            auto cfg = base;
            cfg.Lx = L;
            cfg.inletType = InletType::Uniform;
            auto result = solveCase(cfg);
            outletResults.push_back(result);
            std::cout << "L = " << L
                      << ", iterations = " << result.iterations
                      << ", max y+ = " << result.maxYPlus << std::endl;
        }
        Plotter::plotOutletStudy(lengths, outletResults);*/

        std::cout << "=== Reynolds sweep ===" << std::endl;
        std::vector<double> reValues = {40.0, 80.0, 200, 500, 1000};
        std::vector<double> yPlusValues;
        bool thresholdPrinted = false;

        for (double re : reValues) {
            auto cfg = base;
            cfg.inletType = InletType::Parabolic;
            cfg.Lx = 12.0;
            cfg.reynolds = re;
            cfg.maxIterations = 1200;
            auto result = solveCase(cfg);
            yPlusValues.push_back(result.maxYPlus);
            std::cout << "Re = " << re
                      << ", max y+ = " << result.maxYPlus
                      << ", iterations = " << result.iterations
                      << std::endl;
            if (!thresholdPrinted && result.maxYPlus > 11.63) {
                std::cout << "First case with y+ > 11.63: Re = " << re << std::endl;
                thresholdPrinted = true;
            }
        }
        Plotter::plotReStudy(reValues, yPlusValues);
    }
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_CHANNELSTUDY_H