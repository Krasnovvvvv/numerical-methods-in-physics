#ifndef NUMERICAL_METHODS_IN_PHYSICS_SLAE_EXACTMETHOD_TASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_SLAE_EXACTMETHOD_TASK_H

#pragma once
#include "LabTask.h"
#include "RandomSLAEGenerator.h"
#include <iostream>
#include <vector>

class SLAE_ExactMethod_Task : public LabTask {
public:
    using LabTask::LabTask; // Наследуем конструктор, чтобы ExpectedCurveFunc пробрасывался!

    void run(const std::vector<size_t>& sizes) override {
        std::vector<double> xpoints, ytimes;
        for (auto n : sizes) {
            auto* slaeGen = dynamic_cast<RandomSLAEGenerator*>(&generator);
            if (!slaeGen)
                throw std::runtime_error("Generator must be RandomSLAEGenerator!");

            Eigen::MatrixXd A = slaeGen->generateMatrix(n,true);
            Eigen::VectorXd x_exact = slaeGen->exactSolution(n);
            Eigen::VectorXd b = A * x_exact;

            Timer t;
            auto result = solver.solve(A, b);
            double elapsed = static_cast<double>(t.elapsed());
            xpoints.push_back(static_cast<double>(n));
            ytimes.push_back(elapsed);

            double relerror = (x_exact - result.solution).norm() / x_exact.norm();
            std::cout << n << " | " << relerror << " | " << elapsed << "\n";
        }
        if (sizes.size() > 1)
            plotter.plot(xpoints, ytimes, "Actual", "n", "Time, ms", false, expected);
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_SLAE_EXACTMETHOD_TASK_H
