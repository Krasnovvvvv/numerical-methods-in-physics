#include <iostream>
#include "Labs/Special/Lab2/Base/DichotomousNoiseGenerator.h"
#include "Labs/Special/Lab2/Base/DichotomousRatchetPotential.h"
#include "Labs/Special/Lab2/Tasks/DichotomousRatchetTasks.h"
#include "Labs/Special/Lab2/DiffusionSolvers/LangevinSolver.h"
#include "Labs/Special/Lab2/Base/DiffusionVisualizer.h"

using namespace special;

int main() {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "FULL LABORATORY WORK: DICHOTOMIC RATCHET\n";
    std::cout << "Parts 1 & 2: Noise Generation and Ratchet Dynamics\n";
    std::cout << std::string(80, '=') << "\n";

    // ========================================================================
    // Инициализация
    // ========================================================================

    LangevinSolver solver(42);
    Plotter plotter;

    // Параметры для ЧАСТИ 1 (генерирование шума)
    DiffusionParameters params_noise;
    params_noise.n_particles = 1;
    params_noise.n_steps = 1000;
    params_noise.dt = 0.01;

    // Параметры для ЧАСТИ 2 (рачет-динамика)
    DiffusionParameters params_ratchet;
    params_ratchet.n_particles = 1000;
    params_ratchet.n_steps = 50000;
    params_ratchet.dt = 0.001;
    params_ratchet.diffusion_coeff = 0.1;
    params_ratchet.kB_T = 1.0;
    params_ratchet.x0 = 0.0;
    params_ratchet.x_min = 0.0;
    params_ratchet.x_max = 1.0;
    params_ratchet.use_periodic_bc = true;
    params_ratchet.overdamped = true;

    // ========================================================================
    // ЧАСТЬ 1: ЗАДАНИЯ ПО ГЕНЕРИРОВАНИЮ ДИХОТОМНОГО ШУМА
    // ========================================================================

    std::cout << "\n\n" << std::string(80, '#') << "\n";
    std::cout << "PART 1: DICHOTOMIC NOISE GENERATION\n";
    std::cout << std::string(80, '#') << "\n";

    // Задание 1: Генерирование и проверка свойств шума
    {
        NoiseGenerationTask task1(solver, &plotter, 1,
                                 "Noise Generation");
        task1.run(params_noise);
    }

    // Задание 3: Проверка алгоритма (вычисление моментов)
    {
        NoiseVerificationTask task3(solver, &plotter, 3,
                                   "Noise Verification");
        task3.run(params_noise);
    }

    // ========================================================================
    // ЧАСТЬ 2: ЗАДАНИЯ ПО РАЧЕТ-ПОТЕНЦИАЛАМ
    // ========================================================================

    std::cout << "\n\n" << std::string(80, '#') << "\n";
    std::cout << "PART 2: DICHOTOMIC RATCHET DYNAMICS\n";
    std::cout << std::string(80, '#') << "\n";

    // ЗАДАНИЕ 1A: On-Off рачет с симметричным шумом
    {
        std::cout << "\n>>> Starting Task 1A: On-Off Ratchet (Symmetric Noise)\n";
        OnOffRatchetSymmetric task1a(solver, &plotter, 10, "Task 1A",
                                     1.0,      // V0
                                     1.0,      // L
                                     0.05);    // tau_c
        task1a.run(params_ratchet);
    }

    // ЗАДАНИЕ 2B: Флиппинг рачет с асимметричным шумом
    {
        std::cout << "\n>>> Starting Task 2B: Flipping Ratchet (Asymmetric Noise)\n";
        FlippingRatchetAsymmetric task2b(solver, &plotter, 20, "Task 2B",
                                         1.0,      // V0
                                         1.0,      // L
                                         0.05,     // tau_c
                                         2.0);     // gamma_a / gamma_b
        task2b.run(params_ratchet);
    }

    // ЗАДАНИЕ 2V*: Зависимость скорости от τ_C
    {
        std::cout << "\n>>> Starting Task 2V*: Velocity vs tau_c (Bell-shaped Curve)\n";
        RatchetVelocityVsTauC task2v(solver, &plotter, 30, "Task 2V*",
                                     1.0,      // V0
                                     1.0,      // L
                                     2.0);     // gamma_a / gamma_b
        task2v.run(params_ratchet);
    }

    // ========================================================================
    // Завершение
    // ========================================================================

    std::cout << "\n\n" << std::string(80, '=') << "\n";
    std::cout << "ALL TASKS COMPLETED SUCCESSFULLY\n";
    std::cout << std::string(80, '=') << "\n\n";

    std::cout << "Generated files:\n";
    std::cout << "  - noise_symmetric_tau*.txt     (noise realizations)\n";
    std::cout << "  - noise_asymmetric_tau*.txt    (noise realizations)\n";
    std::cout << "  - ratchet_velocity_vs_tau_c.txt (v(τ_C) dependence)\n";
    std::cout << "\nGenerated plots (via Plotter):\n";
    std::cout << "  - Task 1A: Trajectories, mean position, statistics\n";
    std::cout << "  - Task 2B: Trajectories, mean position, statistics\n";
    std::cout << "  - Task 2V*: Data file for external plotting\n";

    return 0;
}