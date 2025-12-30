#include <iostream>
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"
#include "Helpers/Plotter.h"

int main() {
    std::size_t N    = 500'000;
    std::size_t burn = 10'000;
    double      dt   = 0.01;

    double a = 6.0;
    double b = -4.0;

    std::vector<std::vector<double>> xs, ys;
    std::vector<std::string> labels;

    std::vector<double> tau_cs = {0.5, 2.0, 5.0};

    for (double tau_c : tau_cs) {

        double gamma_b = 1.0 / (2.0 * tau_c);
        double gamma_a = -(a / b) * gamma_b;
        double Gamma   = gamma_a + gamma_b;
        double tau_eff = 1.0 / Gamma;

        std::cout << "tau_c target=" << tau_c
                  << "  actual=" << tau_eff << "\n";

        DichotomicNoise gen(a, b, gamma_a, gamma_b, dt);
        auto prof = gen.generate(N, burn, true);

        // считаем корреляцию до некоторого максимального t'
        double t_max = 40.0;
        std::size_t max_lag = static_cast<std::size_t>(t_max / dt);

        auto C = gen.autocorr_normalized(prof.s, max_lag);

        std::vector<double> tprime(max_lag + 1);
        for (std::size_t k = 0; k <= max_lag; ++k)
            tprime[k] = k * dt;

        xs.push_back(tprime);
        ys.push_back(C);

        labels.push_back("tau_c = " + std::to_string(tau_c));
    }

    Plotter plotter;
    plotter.plot(xs, ys,
                 labels,
                 "t'",
                 "<sigma(t+t') sigma(t)> / <sigma^2>");

}