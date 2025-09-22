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

    // Одиночный и многократный запуск
    void run(const std::vector<size_t>& sizes) override {
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
        if(sizes.size()!=1) {
            plotter->plot(x_points, y_times, "Прогонка", [&y_times](const std::vector<double> &x) {
                // y — экспериментальные значения, должны быть видимы в данной области
                double sum_x2 = 0, sum_xy = 0;
                for (size_t i = 0; i < x.size(); ++i) {
                    sum_x2 += x[i] * x[i];
                    sum_xy += x[i] * y_times[i];
                }
                double alpha = (sum_x2 != 0) ? (sum_xy / sum_x2) : 0.0;

                std::vector<double> y_expected(x.size());
                for (size_t i = 0; i < x.size(); ++i)
                    y_expected[i] = alpha * x[i];
                return y_expected;
            });
        }
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_SLAETASK_H
