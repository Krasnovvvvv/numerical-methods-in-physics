#include "Labs/Special/Additional/NoiseTasks/NoiseMomentTask.h"
#include <iostream>
#include <memory>

int main() {
    double a = 1.0;
    double dt = 0.001;
    double tau = 0.075;
    std::size_t burn_in = 5000;
    std::size_t steps = 5e7;
    std::size_t n_trajectories = 1e2;

    std::unique_ptr<INoiseTask> task = std::make_unique<NoiseMomentTask>(
        a, dt, tau, burn_in, steps, n_trajectories
    );

    std::cout << "Running task: " << task->name() << "\n";

    task -> run();
}