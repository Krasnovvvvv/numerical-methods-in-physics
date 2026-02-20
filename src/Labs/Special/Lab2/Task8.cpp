#include <iostream>
#include <iomanip>
#include "Labs/Special/Additional/RatchetTasks/RatchetTasks.h"

#ifdef _WIN32
#include <windows.h>
#endif

int main() {

#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif

    std::cout << std::setprecision(6) << std::fixed;

    // ==================== ПАРАМЕТРЫ ====================
    double L  = 1.0;         // Период области
    double dt = 0.001;       // Шаг интегрирования
    std::size_t N = 2'000'000;   // Число шагов интегрирования

    /*RatchetTasks::taskA(L, dt, N);

    RatchetTasks::taskB(L, dt, N);

    RatchetTasks::taskC(L, dt, N);

    RatchetTasks::taskG(L, dt, N);

    RatchetTasks::taskD(L, dt, N);*/

    RatchetTasks::taskComparisonAnalyticalPaper(dt, N);

}
