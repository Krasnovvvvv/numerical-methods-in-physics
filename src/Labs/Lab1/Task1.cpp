#include "Labs/Lab1/Tasks/SLAE_ExactMethod_Task.h"
#include "Labs/Lab1/SLAEGenerators/RandomSLAEGenerator.h"
#include "Labs/Lab1/SLAESolvers/ThomasSolver.h"
#include "Helpers/AsymptoticFuncs.h"
#include <vector>

int main() {

    Plotter plot;
    RandomSLAEGenerator gen;
    ThomasSolver solver;

    SLAE_ExactMethod_Task lab(gen,
                              solver,
                              plot,
                              SLAE_ExactMethod_Task::ExpectedCurveFunc(linearAsymptotic));

    std::vector<size_t> sizes;
    for (size_t n = 1000; n <= 2000; n += 100)
        sizes.push_back(n);

    std::cout << std::setw(8) << "n"
              << " | " << std::setw(15) << "relerror"
              << " | " << std::setw(10) << "elapsed (ms)" << std::endl;
    std::cout << std::string(39, '-') << std::endl;

    lab.run(sizes);
}

