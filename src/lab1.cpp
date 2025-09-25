#include "../include/SLAETask.h"
#include "../include/RandomSLAEGenerator.h"
#include "../include/ThomasSolver.h"
#include "../include/AsymptoticFuncs.h"
#include <vector>

int main() {
    Plotter plot;
    RandomSLAEGenerator gen;
    ThomasSolver solver;

    SLAETask lab(gen, solver, plot, SLAETask::ExpectedCurveFunc(linearAsymptotic));
    std::vector<size_t> sizes;
    for (size_t n = 1000; n <= 4000; n += 100)
        sizes.push_back(n);
    lab.run(sizes);

}
