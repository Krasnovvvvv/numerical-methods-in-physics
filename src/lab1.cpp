#include "../include/SLAETask.h"
#include "../include/RandomSLAEGenerator.h"
#include "../include/ThomasSolver.h"
#include <vector>

int main() {

    Plotter plot;
    RandomSLAEGenerator gen;
    ThomasSolver solver;
    SLAETask lab(&gen, &solver, &plot);

    // Эксперимент для графика:
    std::vector<size_t> sizes;
    for (size_t n = 1000; n <= 3000; n += 100) sizes.push_back(n);
    lab.run(sizes);

    return 0;
}
