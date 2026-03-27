#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

#include "Labs/Special/Additional/Base/RatchetParams.h"
#include "Labs/Special/Additional/RatchetTasks/HighTempRatchetTask.h"
#include "Helpers/Plotter.h"

int main() {
    try {
        std::cout << std::setprecision(12);

        RatchetParams params;
        params.V1 = 0.2;
        params.V2 = 0.1;
        params.alpha = -1.0 / 3.0;
        params.epsilon = 0.075;
        params.a = 1.0;
        params.dt = 0.001;

        const std::size_t n_particles = 30000;
        const std::size_t total_time = 500;

        Plotter plotter;

        const unsigned int first_run_seed = 18u;
        const unsigned int n_runs = 2u;

        for (unsigned int i = 0; i < n_runs; ++i) {
            const unsigned int run_seed = first_run_seed + i;

            std::cout << "\n========================================\n";
            std::cout << "Run seed: " << run_seed << "\n";
            std::cout << "========================================\n";

            HighTempRatchetTask task(
                plotter,
                params,
                n_particles,
                total_time,
                false,
                8000,
                run_seed
            );

            task.run();
        }

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << '\n';
        return 2;
    }
}
