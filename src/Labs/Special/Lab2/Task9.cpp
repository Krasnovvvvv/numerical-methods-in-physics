#include <iostream>
#include <exception>

#include "Helpers/Plotter.h"
#include "Labs/Special/Additional/Base/RatchetParams.h"
#include "Labs/Special/Additional/RatchetTasks/HighTempRatchetTask.h"

#ifdef _WIN32
#include <windows.h>
#endif

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif

    try {
        Plotter plotter;

        RatchetParams params;
        params.V1 = 0.05;
        params.V2 = 0.01;
        params.alpha = -1.0 / 3.0;
        params.epsilon = 0.077;
        params.a = 1.0;
        params.dt = 1e-3;

        constexpr std::size_t n_particles = 10000;
        constexpr std::size_t total_time  = 500;

        HighTempRatchetTask task(
            plotter,
            params, n_particles, total_time);
        task.run();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
