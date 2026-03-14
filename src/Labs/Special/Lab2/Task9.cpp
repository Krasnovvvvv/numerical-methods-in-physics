#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

#include "Labs/Special/Additional/Base/RatchetParams.h"
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"
#include "Labs/Special/Additional/RatchetTasks/HighTempRatchetTask.h"
#include "Helpers/Plotter.h"

namespace {

struct NoiseCheckResult {
    double mean_sigma = 0.0;
    double abs_mean_sigma = 0.0;
    double theoretical_mean = 0.0;
    std::size_t n_trajectories = 0;
    std::size_t burn_in_steps = 0;
    std::size_t measure_steps = 0;
};

NoiseCheckResult check_zero_mean_sigma(double a,
                                       double tau_c,
                                       double dt,
                                       std::size_t n_trajectories,
                                       std::size_t burn_in_steps,
                                       std::size_t measure_steps,
                                       unsigned int base_seed = 600u)
{
    NoiseCheckResult res;
    res.n_trajectories = n_trajectories;
    res.burn_in_steps = burn_in_steps;
    res.measure_steps = measure_steps;

    long double sum = 0.0L;
    std::size_t cnt = 0;

    for (std::size_t tr = 0; tr < n_trajectories; ++tr) {
        DichotomicNoise noise(
            a,
            tau_c,
            dt,
            base_seed + static_cast<unsigned int>(tr)
        );

        if (tr == 0) {
            res.theoretical_mean = noise.mean_theoretical();
        }

        for (std::size_t i = 0; i < burn_in_steps; ++i) {
            noise.step();
        }

        for (std::size_t i = 0; i < measure_steps; ++i) {
            sum += static_cast<long double>(noise.step());
            ++cnt;
        }
    }

    if (cnt > 0) {
        res.mean_sigma = static_cast<double>(sum / static_cast<long double>(cnt));
        res.abs_mean_sigma = std::abs(res.mean_sigma);
    }

    return res;
}

void print_noise_check(const NoiseCheckResult& r, double tolerance) {
    std::cout << "=== Zero-mean sigma check ===\n";
    std::cout << "Trajectories      : " << r.n_trajectories << '\n';
    std::cout << "Burn-in steps     : " << r.burn_in_steps << '\n';
    std::cout << "Measure steps     : " << r.measure_steps << '\n';
    std::cout << "Theoretical mean  : " << r.theoretical_mean << '\n';
    std::cout << "Experimental mean : " << r.mean_sigma << '\n';
    std::cout << "|mean|            : " << r.abs_mean_sigma << '\n';
    std::cout << "Tolerance         : " << tolerance << '\n';
}

} // namespace

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

        const std::size_t n_particles = 10000;
        const std::size_t total_time = 200;

        /*
        const std::size_t noise_trajectories = 1000;
        const std::size_t noise_burn_in_steps = 5000;
        const std::size_t noise_measure_steps = 200000;
        const double zero_mean_tolerance = 5e-3;

        const auto check = check_zero_mean_sigma(
            params.a,
            params.epsilon,
            params.dt,
            noise_trajectories,
            noise_burn_in_steps,
            noise_measure_steps
        );

        print_noise_check(check, zero_mean_tolerance);

        if (check.abs_mean_sigma > zero_mean_tolerance) {
            std::cerr << "\n[ERROR] Zero-mean sigma test FAILED.\n";
            std::cerr << "Try increasing noise_measure_steps and/or noise_trajectories,\n";
            std::cerr << "or inspect the dichotomic noise generator.\n";
            return 1;
        }

        std::cout << "[OK] Zero-mean sigma test passed.\n\n";*/

        Plotter plotter;

        HighTempRatchetTask task(
            plotter,
            params,
            n_particles,
            total_time,
            true,
            8000
        );

        task.run();
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << '\n';
        return 2;
    }
}
