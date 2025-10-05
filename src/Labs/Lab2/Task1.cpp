#include "Labs/Lab2/RootFinders/ParallelRootFinder.h"
#include "Labs/Lab2/Tasks/RootFindTask.h"
#include <cmath>

int main() {

    auto func = [](double x) {
        return 2 * std::log(x) + std::sin(std::log(x)) - std::cos(std::log(x));
    };
    auto isInDomain = [](double x) { return x > 0; };

    double a = 0.1, b = 5;
    double step = 1e-4;
    double tolerance = 1e-20;

    ParallelRootFinder finder;

    RootFindTask task(finder);

    task.run(func, isInDomain, a, b, step, tolerance);

    return 0;
}
