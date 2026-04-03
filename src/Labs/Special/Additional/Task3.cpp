#include "Labs/Special/Additional/NoiseTasks/NoiseCorrelationTask.h"
#include "Helpers/Plotter.h"
#include <iostream>
#include <memory>

int main() {
    double a = 1.0;
    double dt = 0.001;
    std::vector tau = {0.01, 0.075, 0.1, 0.5, 0.75};
    std::size_t burn_in = 100;
    std::size_t steps = 5e5;
    std::size_t n_trajectories = 50;
    std::size_t max_lag = std::min(
        steps / 10, static_cast<std::size_t>(10 * 0.8 * 0.25 / dt)
    );

    Plotter plotter;

    std::unique_ptr<INoiseTask> task = std::make_unique<NoiseCorrelationTask>(
        a, dt, tau, burn_in, steps, n_trajectories, max_lag, plotter
    );

    std::cout << "Running task: " << task->name() << "\n";

    task -> run();
}