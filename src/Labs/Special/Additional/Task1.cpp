#include "Labs/Special/Additional/NoiseTasks/NoiseProfileTask.h"
#include "Helpers/Plotter.h"
#include <iostream>

int main() {
    double a  = 1.0;
    double dt = 0.001;
    std::size_t N = 10000;

    std::vector tau_list = {0.5, 0.1, 0.05};

    Plotter plotter;

    std::unique_ptr<INoiseTask> task = std::make_unique<NoiseProfileTask>(
        a, a, dt, N, tau_list, plotter
    );

    std::cout << "Running task: " << task->name() << "\n";

    task->run();
}