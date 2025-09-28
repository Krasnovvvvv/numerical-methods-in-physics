#ifndef NUMERICAL_METHODS_IN_PHYSICS_SLAE_ITERMETHOD_TASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_SLAE_ITERMETHOD_TASK_H

#pragma once
#include "Base/LabTask.h"
#include "Labs/Lab1/SLAEGenerators/RandomSLAEGenerator.h"
#include <vector>

class SLAE_IterMethod_Task : public LabTask {
public:
    using LabTask::LabTask;

    void run(const std::vector<size_t>& sizes) override {
        std::vector<double> xpoints, ypoints;
        auto* slaeGen = dynamic_cast<RandomSLAEGenerator*>(&generator);
        if (!slaeGen) throw std::runtime_error("Generator must be RandomSLAEGenerator!");
        size_t n = sizes.front();
        Eigen::MatrixXd A = slaeGen->generateMatrix(n, true);
        Eigen::VectorXd x_exact = slaeGen->exactSolution(n);
        Eigen::VectorXd b = A * x_exact;

        SolveResult result = solver.solve(A, b);
        for (const auto& [iter, relres] : result.residuals) {
            xpoints.push_back(static_cast<double>(iter));
            ypoints.push_back(relres);
        }
        plotter.plot(xpoints, ypoints, "Relative residual", "Iteration", "||Ax-b||/||b||", true);
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_SLAE_ITERMETHOD_TASK_H
