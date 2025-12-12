#ifndef NUMERICAL_METHODS_IN_PHYSICS_DIFFUSION_PARAMETERS_H
#define NUMERICAL_METHODS_IN_PHYSICS_DIFFUSION_PARAMETERS_H

#pragma once

namespace special {

/**
 * @struct DiffusionParameters
 * @brief Параметры для моделирования диффузии и стохастической динамики
 */
struct DiffusionParameters {
    // ===== ОСНОВНЫЕ ПАРАМЕТРЫ =====

    /// Коэффициент диффузии D
    double diffusion_coeff = 0.1;

    /// Тепловая энергия k_B * T
    double kB_T = 1.0;

    /// Температура (если нужна)
    double T = 1.0;

    /// Постоянная Больцмана (часто = 1 в безразмерных единицах)
    double kB = 1.0;

    /// Постоянная сила
    double constant_force = 0.0;

    // ===== ПАРАМЕТРЫ СЕТКИ И ВРЕМЕНИ =====

    /// Временной шаг dt
    double dt = 0.001;

    /// Количество временных шагов
    int n_steps = 1000;

    /// Период (для периодических граничных условий)
    double L = 1.0;

    /// Левая граница
    double x_min = 0.0;

    /// Правая граница
    double x_max = 1.0;

    // ===== ПАРАМЕТРЫ ЧАСТИЦ =====

    /// Количество частиц
    int n_particles = 100;

    /// Начальное положение частиц
    double x0 = 0.0;

    /// Масса частицы (для инерциального режима)
    double mass = 1.0;

    /// Коэффициент трения γ (для инерциального режима Ланжевена)
    double friction = 1.0;

    // ===== РЕЖИМЫ РЕШЕНИЯ =====

    /// Использовать сверхдемпфированный режим (истина) или полное Ланжевено (ложь)
    bool overdamped = true;

    /// Использовать периодические граничные условия
    bool use_periodic_bc = false;

    // ===== ПАРАМЕТРЫ ВИЗУАЛИЗАЦИИ =====

    /// Количество шагов между кадрами анимации
    int animation_skip_steps = 1;

    /// Максимальное число частиц для отображения в анимации
    int max_particles_to_display = 50;

    // ===== ПОТЕНЦИАЛ =====

    /// Жёсткость гармонического потенциала (U = k*x²/2)
    double potential_stiffness = 1.0;

    /// Амплитуда асимметричного потенциала (рачета)
    double ratchet_amplitude = 0.5;

    /// Волновое число рачета (2π/λ)
    double ratchet_wavenumber = 2.0 * 3.14159265359;

    /// Частота переключения рачета
    double ratchet_frequency = 1.0;

    /// Скорость путешествующей волны
    double traveling_wave_speed = 1.0;
};

} // namespace special

#endif // NUMERICAL_METHODS_IN_PHYSICS_DIFFUSION_PARAMETERS_H