#ifndef NUMERICAL_METHODS_IN_PHYSICS_SLAETASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_SLAETASK_H

#pragma once
#include "LabTask.h"
#include "RandomSLAEGenerator.h"
#include <iostream>
#include <vector>

class SLAETask : public LabTask {
public:
    using LabTask::LabTask;

    // Одиночный запуск
    void run(size_t n) override {
        auto* slaeGen = dynamic_cast<RandomSLAEGenerator*>(generator);
        if (!slaeGen) throw std::runtime_error("Generator must be RandomSLAEGenerator!");

        auto A = slaeGen->generateMatrix(n);
        auto b = slaeGen->generateVector(n);

        Timer<> t;
        auto result = solver->solve(A, b);
        double elapsed = static_cast<double>(t.elapsed());

        auto exact = slaeGen->exactSolution(n);
        double relative_error = (exact - result.solution).norm() / exact.norm();

        std::cout << "Ошибка решения: " << relative_error << std::endl;
        std::cout << "Время решения: " << elapsed << " мс\n";
    }

    void run_experiment(const std::vector<size_t>& sizes) {
        std::vector<double> x_points, y_times;
        for (auto n : sizes) {
            auto* slaeGen = dynamic_cast<RandomSLAEGenerator*>(generator);
            if (!slaeGen) throw std::runtime_error("Generator must be RandomSLAEGenerator!");
            auto A = slaeGen->generateMatrix(n);
            auto b = slaeGen->generateVector(n);

            Timer<> t;
            auto result = solver->solve(A, b);
            double elapsed = static_cast<double>(t.elapsed());

            x_points.push_back(static_cast<double>(n));
            y_times.push_back(elapsed);

            auto exact = slaeGen->exactSolution(n);
            double rel_error = (exact - result.solution).norm() / exact.norm();
            std::cout << "n=" << n << ", Ошибка: " << rel_error << ", Время: " << elapsed << " мс\n";
        }
        plotter->plot_curve(x_points, y_times, "Прогонка");
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_SLAETASK_H
