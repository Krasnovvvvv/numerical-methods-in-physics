#include "Labs/Lab1/SLAESolvers/ThomasSolver.h"
#include "Helpers/Plotter.h"
#include "Labs/Lab6/Tasks/TaskLaserRod.h"

int main() {
    Plotter plotter;
    ThomasSolver solver;

    TaskLaserRod::PhysParams phys {
        .rho = 8000.0,
        .c   = 500.0,
        .k   = 50.0,
        .l   = 0.01,
        .tp  = 0.01,
        .I0  = 1.0e5,
        .T0  = 300.0
    };

    double tMax = 0.05;    // считать до 0.05 c (5 t_p)
    double h_x  = 0.0001;  // шаг по x, 0.1 мм
    double tau  = 0.0001;  // шаг по времени, 0.1 мс

    TaskLaserRod task(solver, &plotter,
                      TaskLaserRod::SchemeType::CrankNicolson);
    task.run(phys, tMax, h_x, tau);

    task.plotTemperatureAtLeft();   // T(0,t)

    return 0;
}
