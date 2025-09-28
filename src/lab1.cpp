#include "../include/SLAE_ExactMethod_Task.h"
#include "../include/RandomSLAEGenerator.h"
#include "../include/ThomasSolver.h"
#include "../include/AsymptoticFuncs.h"
#include <vector>

int main() {
    Plotter plot;
    RandomSLAEGenerator gen;
    ThomasSolver solver;

    SLAE_ExactMethod_Task lab(gen, solver, plot, SLAE_ExactMethod_Task::ExpectedCurveFunc(linearAsymptotic));
    std::vector<size_t> sizes;
    for (size_t n = 1000; n <= 4000; n += 100)
        sizes.push_back(n);
    lab.run(sizes);

}
