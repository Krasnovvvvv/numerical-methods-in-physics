// Lab7_DenseWarmup.cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <functional>
#include <iostream>

#include "Helpers/Plotter.h""
#include "Labs/Lab7/Tasks/TaskWaveExplicit.h"

double exactSolution_piecewiseVelocity(double x, double t,
                                       const TaskWaveBase::PhysParams& phys,
                                       int Nmax = 200)
{
    double c   = phys.c;
    double L   = phys.L;
    double x0  = phys.x0;
    double d   = phys.delta;
    double v0  = phys.v0;

    double a = std::max(0.0, x0 - d);
    double b = std::min(L,  x0 + d);

    if (a >= b)
        return 0.0;

    double sum = 0.0;
    for (int n = 1; n <= Nmax; ++n) {
        double kn = n * M_PI / L;

        // интеграл \int_a^b v0 sin(kn ξ) dξ = v0/kn (cos(kn a) - cos(kn b))
        double In = v0 / kn * (std::cos(kn * a) - std::cos(kn * b));

        // коэффициент B_n = 2/(nπ c)*In = 2/(c L kn) * In
        double Bn = 2.0 / (c * L * kn) * In;

        double sin_kx = std::sin(kn * x);
        double sin_wt = std::sin(kn * c * t);

        sum += Bn * sin_kx * sin_wt;
    }
    return sum;
}

int main() {
    Plotter plotter;

    // физические параметры струны
    TaskWaveBase::PhysParams phys;
    phys.c     = 100.0;   // м/с
    phys.L     = 1.0;     // м
    phys.x0    = 0.3;     // м
    phys.delta = 0.05;    // м
    phys.v0    = 1.0;     // м/с

    double tMax = 0.05;     // с
    double h_x  = 0.01;     // м

    // шаг по времени по условию Куранта (стабильный)
    double tau_stable = 0.5 * h_x / phys.c; // c*tau <= h_x

    // 1. Численное решение явной схемой "крест"
    TaskWaveExplicit taskStable(&plotter);
    taskStable.run(phys, tMax, h_x, tau_stable);

    // 2. u(x0,t) для наглядности
    taskStable.plotAtPoint(phys.x0);

    // 3. Сравнение с точным решением: считаем максимум ошибки
    const auto& x = taskStable.getX();
    const auto& t = taskStable.getTime();
    const auto& u = taskStable.getU();

    double maxErr = 0.0;
    for (int s = 0; s < t.size(); ++s) {
        for (int i = 0; i < x.size(); ++i) {
            double u_ex = exactSolution_piecewiseVelocity(x[i], t[s], phys);
            double diff = std::abs(u(s, i) - u_ex);
            if (diff > maxErr) maxErr = diff;
        }
    }
    std::cout << "Max error vs exact (explicit scheme, CFL satisfied): "
              << maxErr << "\n";

    // 4. Демонстрация нарушения условия Куранта
    double tau_unstable = 2.0 * h_x / phys.c; // c*tau > h_x

    TaskWaveExplicit taskUnstable(&plotter);
    taskUnstable.run(phys, tMax, h_x, tau_unstable);

    taskUnstable.plotAtPoint(phys.x0);

    std::cout << "Stable tau = "   << tau_stable
              << ", unstable tau = " << tau_unstable << "\n";

    return 0;
}
