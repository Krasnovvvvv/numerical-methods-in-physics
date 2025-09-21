#include "../include/SLAETask.h"
#include "../include/RandomSLAEGenerator.h"
#include "../include/ThomasSolver.h"
#include <vector>

int main() {
    // Ожидаемая кривая
    auto expected_O_n = [](const std::vector<double>& x) {
        double alpha = 0.01; // подобрать по эксперименту
        std::vector<double> y(x.size());
        for (size_t i = 0; i < x.size(); ++i) y[i] = alpha * x[i];
        return y;
    };

    Plotter plot(std::make_optional(expected_O_n));
    RandomSLAEGenerator gen;
    ThomasSolver solver;
    SLAETask lab(&gen, &solver, &plot);

    // Эксперимент для графика:
    std::vector<size_t> sizes;
    for (size_t n = 1000; n <= 3000; n += 100) sizes.push_back(n);
    lab.run_experiment(sizes);

    return 0;
}
