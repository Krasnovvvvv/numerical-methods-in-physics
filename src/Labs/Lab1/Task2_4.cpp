#include "Labs/Lab1/Tasks/SLAE_IterMethod_Task.h"
#include "Labs/Lab1/SLAEGenerators/RandomSLAEGenerator.h"
#include "Labs/Lab1/SLAESolvers/GradientDescentSolver.h"

int main() {

    Plotter plot;
    RandomSLAEGenerator gen(-1000,1000);
    GradientDescentSolver solver(1000,1e-5);

    SLAE_IterMethod_Task task(gen,solver,plot);
    task.run({1000});
}




