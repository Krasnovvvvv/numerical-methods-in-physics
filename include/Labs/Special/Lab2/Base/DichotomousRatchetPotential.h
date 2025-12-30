#ifndef DICHOTOMIC_RATCHET_POTENTIAL_H
#define DICHOTOMIC_RATCHET_POTENTIAL_H

#pragma once

#include "Potential1D.h"
#include "DichotomousNoiseGenerator.h"
#include <cmath>
#include <vector>
#include <memory>
#include <iostream>

namespace special {

/**
 * @class DichotomousRatchetPotential
 * @brief Асимметричный потенциал (рачет) с дихотомным модулированием
 *
 * - On-Off рачет: U(x,t) = V₀·V(x)·[(1+σ(t))/2]
 * - Flipping рачет: U(x,t) = V₀·V(x)·[ε + σ(t)]
 *
 * Базовый бигармонический потенциал (уравнение E1):
 * V(x) = (1/2π)·[sin(2πx/L) + A₂·sin(4πx/L)]
 */
class DichotomousRatchetPotential : public Potential1D {
public:
    /// Тип модуляции потенциала
    enum class ModulationType {
        ON_OFF,        ///< f(t) = (1 + σ)/2, σ ∈ {-1, 1} → f ∈ {0, 1}
        FLIPPING       ///< f(t) = ε + σ, σ ∈ {-1, 1} → f ∈ {ε-1, ε+1}
    };

    /**
     * @brief Конструктор
     * @param V0 амплитуда потенциала
     * @param L период потенциала (обычно 1.0)
     * @param A2 коэффициент второй гармоники (обычно 0.25 для асимметрии)
     * @param gamma_a частота переходов a→b дихотомного шума
     * @param gamma_b частота переходов b→a дихотомного шума
     * @param modulation_type тип модуляции (ON_OFF или FLIPPING)
     * @param epsilon параметр для FLIPPING режима (обычно 0.0)
     * @param seed для генератора случайных чисел
     */
    DichotomousRatchetPotential(
        double V0,
        double L,
        double A2,
        double gamma_a = 1.0,
        double gamma_b = 1.0,
        ModulationType modulation_type = ModulationType::ON_OFF,
        double epsilon = 0.0,
        unsigned int seed = 42)
        : V0_(V0), L_(L), A2_(A2),
          gamma_a_(gamma_a), gamma_b_(gamma_b),
          modulation_type_(modulation_type),
          epsilon_(epsilon),
          current_sigma_(1.0),  // начальное значение
          current_time_(-1.0),
          generator_(std::make_shared<DichotomousNoiseGenerator>(
              1.0, -1.0, gamma_a, gamma_b, seed))
    {
    }

    /**
     * @brief Установить текущее значение σ(t) из внешнего генератора
     * Вызывать перед вычислением U и F в каждый момент времени
     * @param sigma текущее значение дихотомного шума
     * @param t текущее время
     */
    void set_sigma(double sigma, double t) {
        current_sigma_ = sigma;
        current_time_ = t;
    }

    /**
     * @brief Получить текущее значение σ(t)
     */
    double get_sigma() const {
        return current_sigma_;
    }

    /**
     * @brief Базовый потенциал V(x) (без модуляции, без множителя V₀)
     *
     * V(x) = (1/2π)·[sin(2πx/L) + A₂·sin(4πx/L)]
     *
     * уравнение E1 из лекции
     */
    double V_base(double x) const {
        const double arg1 = 2.0 * M_PI * x / L_;
        const double arg2 = 4.0 * M_PI * x / L_;
        return (std::sin(arg1) + A2_ * std::sin(arg2)) / (2.0 * M_PI);
    }

    /**
     * @brief Производная базового потенциала V'(x)
     *
     * dV/dx = cos(2πx/L) + (A₂/2)·cos(4πx/L)
     */
    double V_base_prime(double x) const {
        const double arg1 = 2.0 * M_PI * x / L_;
        const double arg2 = 4.0 * M_PI * x / L_;
        return (std::cos(arg1) + A2_ * std::cos(arg2));
    }

    /**
     * @brief Функция модуляции f(t)
     *
     * ON_OFF:  f(t) = (1 + σ(t))/2  → f ∈ {0, 1}
     * FLIPPING: f(t) = ε + σ(t)     → f ∈ {ε-1, ε+1}
     */
    double modulation_f() const {
        switch (modulation_type_) {
            case ModulationType::ON_OFF:
                // f(t) = (1 + σ(t))/2
                // σ ∈ {-1, 1} ⇒ f ∈ {0, 1}
                return 0.5 * (1.0 + current_sigma_);

            case ModulationType::FLIPPING:
                // f(t) = ε + σ(t)
                // σ ∈ {-1, 1} ⇒ f ∈ {ε-1, ε+1}
                return epsilon_ + current_sigma_;

            default:
                return 1.0;
        }
    }

    /**
     * @brief Потенциальная энергия U(x,t)
     *
     * U(x,t) = V₀·V(x)·f(t)
     *
     * где V(x) — базовый потенциал
     * f(t) — временная модуляция через дихотомный шум σ(t)
     */
    double U(double x, double t = 0.0) const override {
        return V0_ * V_base(x) * modulation_f();
    }

    /**
     * @brief Сила F(x,t) = -dU/dx
     *
     * F = -dU/dx = -V₀·V'(x)·f(t)
     */
    double F(double x, double t = 0.0) const override {
        // F = -dU/dx = -V₀·V'(x)·f(t)
        return -V0_ * V_base_prime(x) * modulation_f();
    }

    /**
     * @brief Название потенциала
     */
    std::string name() const override {
        std::string type_str = (modulation_type_ == ModulationType::ON_OFF)
            ? "OnOff" : "Flipping";
        return "DichotomousRatchet[" + type_str + "]";
    }

    /**
     * @brief Получить доступ к генератору дихотомного шума
     */
    std::shared_ptr<DichotomousNoiseGenerator> generator() {
        return generator_;
    }

    /**
     * @brief Вывести параметры потенциала
     */
    void print_parameters() const {
        std::cout << "\n=== DichotomousRatchetPotential Parameters ===\n";
        std::cout << "V₀ = " << V0_ << "\n";
        std::cout << "L = " << L_ << "\n";
        std::cout << "A₂ = " << A2_ << "\n";
        std::cout << "γ_a = " << gamma_a_ << ", γ_b = " << gamma_b_ << "\n";
        std::cout << "Γ = " << (gamma_a_ + gamma_b_) << "\n";
        std::cout << "τ_C = " << (1.0 / (gamma_a_ + gamma_b_)) << "\n";
        std::cout << "Modulation: "
                  << (modulation_type_ == ModulationType::ON_OFF ? "ON_OFF" : "FLIPPING")
                  << "\n";
        if (modulation_type_ == ModulationType::FLIPPING) {
            std::cout << "ε = " << epsilon_ << "\n";
        }
        std::cout << "=============================================\n\n";
    }

private:
    double V0_;           // амплитуда потенциала
    double L_;            // период потенциала
    double A2_;           // коэффициент асимметрии
    double gamma_a_, gamma_b_;  // частоты переходов дихотомного шума
    ModulationType modulation_type_;
    double epsilon_;      // смещение для FLIPPING режима
    mutable double current_sigma_;  // текущее значение σ(t)
    mutable double current_time_;   // текущее время
    std::shared_ptr<DichotomousNoiseGenerator> generator_;
};

/**
 * @class SymmetricRatchetPotential
 * @brief Контрольный потенциал: пилообразный без асимметрии
 *
 * Используется для проверки: при симметричном потенциале
 * и симметричном дихотомном шуме должно быть v = 0
 */
class SymmetricRatchetPotential : public Potential1D {
public:
    SymmetricRatchetPotential(double V0 = 1.0, double L = 1.0)
        : V0_(V0), L_(L)
    {
    }

    /**
     * @brief Чистый синусоидальный потенциал (симметричный)
     * U(x) = V₀·sin²(πx/L)
     */
    double U(double x, double t = 0.0) const override {
        double s = std::sin(M_PI * x / L_);
        return V0_ * s * s;
    }

    /**
     * @brief Производная
     * F = -dU/dx = -V₀·2·sin(πx/L)·cos(πx/L)·(π/L)
     *             = -V₀·sin(2πx/L)·(π/L)
     */
    double F(double x, double t = 0.0) const override {
        return -V0_ * std::sin(2.0 * M_PI * x / L_) * (M_PI / L_);
    }

    std::string name() const override {
        return "SymmetricRatchet";
    }

private:
    double V0_, L_;
};

} // namespace special

#endif // DICHOTOMIC_RATCHET_POTENTIAL_H