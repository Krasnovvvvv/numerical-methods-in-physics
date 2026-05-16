#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "Labs/Special/Lab7/Task/AdvectionDiffusionTask.h"

namespace {
    constexpr double X0 = 0.0;
    constexpr double Y0 = 1.0;
    constexpr double X1 = 1.0;
    constexpr double Y1 = 0.0;
    constexpr size_t LinePoints = 401;

    double exactTask2(double x, double y) {
        if (y <= x) {
            return 0.0;
        }

        const double leftValue = ((y - x) <= 0.5) ? 1.0 : 0.0;
        return leftValue;
    }

    AdvectionDiffusionTaskConfig makeTask2Config(double gamma,
                                                 ConvectionScheme scheme,
                                                 TVDLimiterType limiter,
                                                 const std::string& title,
                                                 bool showHeatmap = false,
                                                 bool showHistory = false) {
        AdvectionDiffusionTaskConfig cfg;

        cfg.title = title;

        cfg.Nx = 160;
        cfg.Ny = 160;
        cfg.Lx = 1.0;
        cfg.Ly = 1.0;

        cfg.rho = 1.0;
        cfg.u = 2.0;
        cfg.v = 2.0;
        cfg.Gamma = gamma;

        cfg.Su = 0.0;
        cfg.Sp = 0.0;

        cfg.left = BoundaryCondition::dirichlet([](double y) {
            return (y <= 0.5) ? 1.0 : 0.0;
        });
        cfg.bottom = BoundaryCondition::dirichlet(0.0);
        cfg.right  = BoundaryCondition::outlet();
        cfg.top    = BoundaryCondition::outlet();

        cfg.initialValue = 0.0;

        cfg.scheme = scheme;
        cfg.limiter = limiter;
        cfg.swebyBeta = 1.5;

        cfg.linearSolverType = EigenSolverType::SparseLU;
        cfg.linearMaxIterations = 5000;
        cfg.linearTolerance = 1e-12;

        cfg.outerMaxIterations = 100;
        cfg.outerTolerance = 1e-10;

        cfg.showHeatmap = showHeatmap;
        cfg.showProfileX = false;
        cfg.showProfileY = false;
        cfg.showHistory = showHistory;

        cfg.profileY = 0.50;
        cfg.profileX = 0.50;

        return cfg;
    }

    std::vector<double> exactOnSegmentTask2(size_t points) {
        std::vector<double> exact(points, 0.0);
        for (size_t k = 0; k < points; ++k) {
            const double t = static_cast<double>(k) / static_cast<double>(points - 1);
            const double x = (1.0 - t) * X0 + t * X1;
            const double y = (1.0 - t) * Y0 + t * Y1;
            exact[k] = exactTask2(x, y);
        }
        return exact;
    }

    bool updateOk(bool ok, const AdvectionDiffusionTaskResult& r) {
        return ok && r.solve.converged;
    }
}

int main() {
    bool ok = true;

    const std::vector<TVDLimiterType> tvdLimiters = {
        TVDLimiterType::VanLeer,
        TVDLimiterType::VanAlbada,
        TVDLimiterType::MinMod,
        TVDLimiterType::Superbee,
        TVDLimiterType::UMIST
    };

    auto mainCfg = makeTask2Config(
        0.0,
        ConvectionScheme::TVD,
        TVDLimiterType::VanLeer,
        "Task 2: c(x,y), Gamma = 0, TVD-VanLeer",
        true,
        true
    );
    auto mainResult = AdvectionDiffusionTask::run(mainCfg);
    ok = updateOk(ok, mainResult);

    auto upwindCfg = makeTask2Config(
        0.0,
        ConvectionScheme::Upwind,
        TVDLimiterType::VanLeer,
        "Task 2: Upwind, Gamma = 0"
    );
    auto upwindResult = AdvectionDiffusionTask::solve(upwindCfg);
    ok = updateOk(ok, upwindResult);

    std::vector<std::vector<double>> curves;
    std::vector<std::string> labels;

    auto [lineX, upwindLine] = CFDPlotter::sampleSegment(
        upwindResult.field, X0, Y0, X1, Y1, LinePoints
    );

    curves.push_back(exactOnSegmentTask2(LinePoints));
    labels.emplace_back("Exact");

    curves.push_back(upwindLine);
    labels.emplace_back("Upwind");

    for (auto lim : tvdLimiters) {
        auto cfg = makeTask2Config(
            0.0,
            ConvectionScheme::TVD,
            lim,
            "Task 2: " + AdvectionDiffusionTask::schemeName(ConvectionScheme::TVD, lim) + ", Gamma = 0"
        );

        auto res = AdvectionDiffusionTask::solve(cfg);
        ok = updateOk(ok, res);

        auto [s, c] = CFDPlotter::sampleSegment(
            res.field, X0, Y0, X1, Y1, LinePoints
        );

        curves.push_back(c);
        labels.push_back(AdvectionDiffusionTask::limiterName(lim));
    }

    CFDPlotter::plotCurves(
        lineX,
        curves,
        labels,
        "Task 2: comparison on line (0,1)-(1,0), Gamma = 0",
        "t on line: (x,y) = (t, 1 - t)",
        "c"
    );

    auto diffCfg1 = makeTask2Config(
        1e-3,
        ConvectionScheme::TVD,
        TVDLimiterType::VanLeer,
        "Task 2: c(x,y), Gamma = 1e-3, TVD-VanLeer",
        true,
        false
    );
    auto diffRes1 = AdvectionDiffusionTask::run(diffCfg1);
    ok = updateOk(ok, diffRes1);

    auto diffCfg2 = makeTask2Config(
        1e-2,
        ConvectionScheme::TVD,
        TVDLimiterType::VanLeer,
        "Task 2: c(x,y), Gamma = 1e-2, TVD-VanLeer",
        true,
        false
    );
    auto diffRes2 = AdvectionDiffusionTask::run(diffCfg2);
    ok = updateOk(ok, diffRes2);

    auto [s0, c0] = CFDPlotter::sampleSegment(mainResult.field, X0, Y0, X1, Y1, LinePoints);
    auto [s1, c1] = CFDPlotter::sampleSegment(diffRes1.field, X0, Y0, X1, Y1, LinePoints);
    auto [s2, c2] = CFDPlotter::sampleSegment(diffRes2.field, X0, Y0, X1, Y1, LinePoints);

    CFDPlotter::plotCurves(
        s0,
        {c0, c1, c2},
        {"Gamma = 0", "Gamma = 1e-3", "Gamma = 1e-2"},
        "Task 2: effect of nonzero diffusion on line (0,1)-(1,0)",
        "t on line: (x,y) = (t, 1 - t)",
        "c"
    );

    std::cout << "\nTask 2 complete.\n";
    return ok ? 0 : 1;
}