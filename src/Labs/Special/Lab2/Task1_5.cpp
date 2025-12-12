#include "Labs/Special/Lab2/Base/DiffusionParameters.h"
#include "Labs/Special/Lab2/Base/Potential1D.h"
#include "Labs/Special/Lab2/Base/BaseDiffusionSolver.h"
#include "Labs/Special/Lab2/DiffusionSolvers/LangevinSolver.h"
#include "Labs/Special/Lab2/Base/StochasticTask.h"
#include "Labs/Special/Lab2/Tasks/StochasticTasks.h"
#include "Helpers/Plotter.h"

#include <iostream>
#include <string>

using namespace special;

int main() {
    // ===== ИНИЦИАЛИЗАЦИЯ =====
    Plotter plotter;  // ← ЗДЕСЬ СОЗДАЁМ PLOTTER!
    LangevinSolver solver(42);

    DiffusionParameters params;
    params.L = 1.0;
    params.dt = 0.001;
    params.n_particles = 1000;
    params.x0 = 0.0;
    params.use_periodic_bc = false;
    params.x_min = -10.0;
    params.x_max = 10.0;
    params.overdamped = true;

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "СТОХАСТИЧЕСКАЯ ДИНАМИКА: ЗАДАЧИ 1-5\n";
    std::cout << std::string(70, '=') << "\n";

    // =========================================================================
    // ЗАДАЧА 1: СВОБОДНАЯ ДИФФУЗИЯ
    // =========================================================================
    {
        DiffusionParameters p1 = params;
        p1.diffusion_coeff = 0.1;
        p1.kB_T = 1.0;
        p1.n_steps = 20000;

        FreeDiffusion task1(solver, &plotter, 1, "Free Diffusion");
        task1.run(p1);
    }

    // =========================================================================
    // ЗАДАЧА 2: АНАЛИЗ МОМЕНТОВ
    // =========================================================================
    {
        DiffusionParameters p2 = params;
        p2.diffusion_coeff = 0.1;
        p2.kB_T = 1.0;
        p2.n_steps = 20;
        p2.n_particles = 1000;

        MomentsAnalysis task2(solver, &plotter, 2, "Moments Analysis");
        task2.run(p2);
    }

    // =========================================================================
    // ЗАДАЧА 3: ДИФФУЗИЯ В ПОТЕНЦИАЛЕ
    // =========================================================================
    {
        DiffusionParameters p3 = params;
        p3.diffusion_coeff = 0.1;
        p3.kB_T = 1.0;
        p3.n_steps = 50000;

        DiffusionInPotential task3(solver, &plotter, 3, "Diffusion in Potential");
        task3.run(p3);
    }

    // =========================================================================
    // ЗАДАЧА 4: МИГАЮЩИЙ РАЧЕТ
    // =========================================================================
    {
        DiffusionParameters p4 = params;
        p4.diffusion_coeff = 0.05;
        p4.kB_T = 1.0;
        p4.n_steps = 50000;
        p4.use_periodic_bc = true;

        RatchetFlashing task4(solver, &plotter, 4, "Flashing Ratchet");
        task4.run(p4);
    }

    // =========================================================================
    // ЗАДАЧА 5: НАКЛОННЫЙ РАЧЕТ
    // =========================================================================
    {
        DiffusionParameters p5 = params;
        p5.diffusion_coeff = 0.05;
        p5.kB_T = 1.0;
        p5.n_steps = 50000;
        p5.use_periodic_bc = true;

        RatchetTilting task5(solver, &plotter, 5, "Tilting Ratchet");
        task5.run(p5);
    }

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "✓ ALL TASKS COMPLETED SUCCESSFULLY!\n";
    std::cout << std::string(70, '=') << "\n\n";

    return 0;
}