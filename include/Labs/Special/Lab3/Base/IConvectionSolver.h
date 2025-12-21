#ifndef NUMERICAL_METHODS_IN_PHYSICS_ICONVECTIONSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_ICONVECTIONSOLVER_H

#pragma once

#include <vector>
#include <string>

/// Интерфейс решателя для 1D уравнения переноса:
/// U_t + u U_x = 0, c = u*dt/dx.
struct IConvectionSolver {
    virtual ~IConvectionSolver() = default;

    /// Один шаг по времени: U^n -> U^{n+1}.
    /// u_cur  — значения на текущем слое,
    /// u_next — новый слой,
    /// c      — конвекционное число u*dt/dx.
    virtual void step(std::vector<double>& u_next,
                      const std::vector<double>& u_cur,
                      double c) const = 0;

    [[nodiscard]] virtual std::string name() const = 0;
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_ICONVECTIONSOLVER_H
