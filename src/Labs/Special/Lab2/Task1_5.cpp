#include "Labs/Special/Lab2/Base/DiffusionParameters.h"
#include "Labs/Special/Lab2/Base/Potential1D.h"
#include "Labs/Special/Lab2/Base/BaseDiffusionSolver.h"
#include "Labs/Special/Lab2/DiffusionSolvers/LangevinSolver.h"
#include "Labs/Special/Lab2/Base/StochasticTask.h"
#include "Labs/Special/Lab2/Tasks/StochasticTasks.h"
#include "Helpers/Plotter.h"

#include <iostream>
#include <string>
#include <cmath>
#include <clocale>
#include <windows.h>

using namespace special;

int main() {
    setlocale(LC_ALL, "");
#ifdef _WIN32
#include <windows.h>
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif

    // ===== ИНИЦИАЛИЗАЦИЯ =====
    Plotter plotter;
    LangevinSolver solver(42);

    DiffusionParameters base;
    base.L              = 1.0;
    base.dt             = 0.001;
    base.n_particles    = 1000;
    base.x0             = 0.0;
    base.use_periodic_bc = false;
    base.x_min          = -10.0;
    base.x_max          = 10.0;
    base.overdamped     = true;
    base.diffusion_coeff = 0.1;
    base.kB_T           = 1.0;

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "СТОХАСТИЧЕСКАЯ ДИНАМИКА: ЗАДАЧИ 1-5\n";
    std::cout << std::string(70, '=') << "\n";

    // =========================================================================
    // ЗАДАЧА 1: СВОБОДНАЯ ДИФФУЗИЯ
    // =========================================================================
    {
        DiffusionParameters p1 = base;
        p1.n_steps = 20000;

        FreeDiffusion task1(solver, &plotter, 1, "Task 1: Free diffusion");
        task1.run(p1);
    }

    // =========================================================================
    // ЗАДАЧА 2: АНАЛИЗ МОМЕНТОВ (N = 20)
    // =========================================================================
    {
        DiffusionParameters p2 = base;
        p2.n_steps     = 20;
        p2.n_particles = 1000;

        MomentsAnalysis task2(solver, &plotter, 2, "Task 2: Moments");
        task2.run(p2);
    }

    // =========================================================================
    // ЗАДАЧА 3: (уже покрыта визуализацией в задаче 1 — распределения с теорией)
    // =========================================================================

    // =========================================================================
    // ЗАДАЧА 4: СТАЦИОНАРНАЯ ОДНОРОДНАЯ СИЛА  F(x,t) = F
    // =========================================================================
    {
        DiffusionParameters p4 = base;
        p4.n_steps        = 20000;
        p4.constant_force = 5;

        ConstantForceDiffusion task4(
            solver, &plotter, 4, "Task 4: constant force");
        task4.run(p4);
    }

    // Параметры периодического потенциала для задач 5A и 5B
    double V0 = 1.0;
    double Lp = 2.0 * M_PI;
    double F  = 0.5;

    // =========================================================================
    // ЗАДАЧА 5A: ДИФФУЗИЯ В СТАЦИОНАРНОМ ПЕРИОДИЧЕСКОМ ПОТЕНЦИАЛЕ
    // U(x) = V0 sin(2π x / L)
    // =========================================================================
    {
        DiffusionParameters p5a = base;
        p5a.n_steps        = 50000;
        p5a.use_periodic_bc = true;

        PeriodicPotentialDiffusion task5a(
            solver, &plotter, 5, "Task 5A: periodic potential", V0, Lp);
        task5a.run(p5a);
    }

    // =========================================================================
    // ЗАДАЧА 5B: ТОТ ЖЕ ПЕРИОДИЧЕСКИЙ ПОТЕНЦИАЛ ПОД ДЕЙСТВИЕМ ПОСТОЯННОЙ СИЛЫ
    // U(x) = V0 sin(2π x / L) - F x
    // =========================================================================
    {
        DiffusionParameters p5b = base;
        p5b.n_steps        = 50000;
        p5b.use_periodic_bc = true;

        PeriodicPotentialWithForce task5b(
            solver, &plotter, 6, "Task 5B: periodic + constant force",
            V0, Lp, F);
        task5b.run(p5b);
    }

    // =========================================================================
    // ЗАДАЧА 5С: Бегущий потенциал
    // U(x,t) = V0 * sin^2[π (x - v t) / L]
    // =========================================================================
    {
        DiffusionParameters p5c = base;
        p5c.diffusion_coeff  = 0.05;
        p5c.kB_T             = 1.0;
        p5c.n_steps          = 30000;
        p5c.use_periodic_bc  = true;
        p5c.traveling_wave_speed = 2.5;

        // Параметры потенциала
        double V0 = 4.0;           // амплитуда барьеров
        double L  = 1.0;           // период по x

        TravelingWaveRatchet task5c(
            solver,
            &plotter,
            6,
            "Traveling-wave ratchet",
            V0, L);

        task5c.run(p5c);
    }

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "✓ ALL TASKS COMPLETED SUCCESSFULLY!\n";
    std::cout << std::string(70, '=') << "\n\n";

    return 0;
}