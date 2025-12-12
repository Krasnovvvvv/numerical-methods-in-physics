#ifndef NUMERICAL_METHODS_IN_PHYSICS_BASE_DIFFUSION_SOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_BASE_DIFFUSION_SOLVER_H

#pragma once

#include "DiffusionParameters.h"
#include "Potential1D.h"
#include <vector>
#include <string>
#include <memory>

namespace special {

/**
 * @struct DiffusionResult
 * @brief Результат решения задачи диффузии
 */
struct DiffusionResult {
    std::vector<double> t;           ///< Моменты времени
    std::vector<double> x;           ///< Положения частиц (для Ланжевена)
    std::vector<std::vector<double>> trajectories;  ///< Траектории частиц [particle][time]
    std::vector<std::vector<double>> rho;           ///< Плотность вероятности [time][x]
    std::vector<double> J;           ///< Поток вероятности
    std::vector<double> mean_x;      ///< Среднее положение <x(t)>
    std::vector<double> mean_x2;     ///< <x²(t)>
    std::vector<double> mean_x3;     ///< <x³(t)> – третий момент (НОВОЕ)
    std::vector<double> mean_x4;     ///< <x⁴(t)> – четвёртый момент (НОВОЕ)
    std::vector<double> displacement_squared;  ///< <(x-x0)²>
    std::vector<double> velocity;    ///< Средняя скорость <v(t)>
    int steps{0};
};

/**
 * @class BaseDiffusionSolver
 * @brief Абстрактный базовый класс для решателей диффузионной динамики
 */
class BaseDiffusionSolver {
public:
    virtual ~BaseDiffusionSolver() = default;

    /**
     * @brief Основной метод решения
     * @param params Параметры задачи
     * @param potential Потенциальное поле
     * @return Результат расчёта
     */
    virtual DiffusionResult solve(const DiffusionParameters& params,
                                  const Potential1D& potential) = 0;

    /**
     * @brief Имя решателя для вывода
     */
    virtual std::string name() const = 0;

    /**
     * @brief Сделать один шаг временной интеграции (для анимации)
     * @param params Параметры
     * @param potential Потенциал
     * @param state Текущее состояние [x_1, x_2, ..., x_n] для n частиц
     * @param time Текущее время
     * @return Новое состояние
     */
    virtual std::vector<double> step(const DiffusionParameters& params,
                                      const Potential1D& potential,
                                      const std::vector<double>& state,
                                      double time) = 0;
};

} // namespace special

#endif // NUMERICAL_METHODS_IN_PHYSICS_BASE_DIFFUSION_SOLVER_H