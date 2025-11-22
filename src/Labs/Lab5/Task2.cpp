#include "Labs/Lab1/SLAESolvers/ThomasSolver.h"
#include "Labs/Lab5/Tasks/HeatRodTask.h"
#include <cmath>

int main() {

    std::vector<int> Nodes = {500};
    std::vector<double> hP_over_kA = {22.0, 77.0};
    double L = 0.1;
    double Tb = 420.0;
    double T_inf = 290.0;

    ThomasSolver solver;
    Plotter plotter;
    HeatRodTask task(solver, &plotter);

    task.run(hP_over_kA, L, Nodes, Tb, T_inf);
}