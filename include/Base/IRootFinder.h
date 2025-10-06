#ifndef NUMERICAL_METHODS_IN_PHYSICS_IROOTFINDER_H
#define NUMERICAL_METHODS_IN_PHYSICS_IROOTFINDER_H

#pragma once
#include <functional>
#include <optional>

// Результат: найденный интервал смены знака
struct RootInterval {
    double left;
    double right;
};

class IRootFinder {
public:
    virtual ~IRootFinder() = default;

    // f — функция, isInDomain — область допустимых значений
    virtual std::optional<RootInterval> findSignChange(
        std::function<double(double)> f,
        std::function<bool(double)> isInDomain, // функция ОДЗ
        double a, double b, double step, double tol
    ) = 0;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_IROOTFINDER_H