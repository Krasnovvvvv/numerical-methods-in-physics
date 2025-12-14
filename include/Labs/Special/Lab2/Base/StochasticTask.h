#ifndef NUMERICAL_METHODS_IN_PHYSICS_STOCHASTIC_TASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_STOCHASTIC_TASK_H

#pragma once

#include "BaseDiffusionSolver.h"
#include "DiffusionParameters.h"
#include "Potential1D.h"
#include "Labs/Special/Lab2/Base/DiffusionVisualizer.h"
#include "Helpers/Plotter.h"
#include <memory>
#include <iostream>
#include <iomanip>

namespace special {

/**
 * @class StochasticTask
 * @brief Базовый класс для задач диффузии с визуализацией
 */
class StochasticTask {
public:
    StochasticTask(BaseDiffusionSolver& solver,
                   Plotter* plotter,
                   unsigned short task_id,
                   const std::string& task_name)
        : solver_(solver),
          plotter_(plotter),
          task_id_(task_id),
          task_name_(task_name),
          visualizer_(plotter) {}

    virtual ~StochasticTask() = default;

    /**
     * @brief Запустить полную задачу
     */
    void run(const DiffusionParameters& params) {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "Task #" << task_id_ << ": " << task_name_ << "\n";
        std::cout << std::string(70, '=') << "\n";

        // 1. Создать потенциал
        auto potential = create_potential(params);

        // 2. Решить задачу
        std::cout << "Solver: " << solver_.name() << "\n";
        std::cout << "Particles: " << params.n_particles << ", Steps: " << params.n_steps << "\n";
        DiffusionResult result = solver_.solve(params, *potential);
        std::cout << "✓ Solved!\n\n";

        // 3. Постобработка
        postprocess(result, params);

        // 4. Визуализация
        std::cout << "Plotting graphs...\n";
        visualize(result, params);
        std::cout << "\n";
    }

    /**
     * @brief Вывести статистику
     */
    void print_statistics(const DiffusionResult& result,
                         const DiffusionParameters& params) {
        if (result.t.empty() || result.mean_x.empty()) return;

        double final_mean_x = result.mean_x.back();
        double final_disp = result.displacement_squared.back();
        double final_time = result.t.back();

        std::cout << "Statistics:\n";
        std::cout << "  <x>: " << std::scientific << final_mean_x << "\n";
        std::cout << "  <(x-x0)²>: " << final_disp << "\n";
        std::cout << "  Time: " << final_time << "\n";

        if (final_time > 1e-6) {
            double estimated_D = final_disp / (2.0 * final_time);
            std::cout << "  Estimated D: " << estimated_D << "\n";
            std::cout << "  True D: " << params.diffusion_coeff << "\n";
        }

        std::cout << std::defaultfloat << "\n";
    }

protected:
    BaseDiffusionSolver& solver_;
    Plotter* plotter_;
    unsigned short task_id_;
    std::string task_name_;
    DiffusionVisualizer visualizer_;

    virtual std::unique_ptr<Potential1D> create_potential(const DiffusionParameters& params) = 0;

    virtual void postprocess(const DiffusionResult& result,
                            const DiffusionParameters& params) {
        print_statistics(result, params);
    }

    virtual void visualize(const DiffusionResult& result,
                          const DiffusionParameters& params) {
        // Переопределить в подклассах
    }
};

} // namespace special

#endif // NUMERICAL_METHODS_IN_PHYSICS_STOCHASTIC_TASK_H