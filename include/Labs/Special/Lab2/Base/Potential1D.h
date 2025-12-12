#ifndef NUMERICAL_METHODS_IN_PHYSICS_POTENTIAL_1D_H
#define NUMERICAL_METHODS_IN_PHYSICS_POTENTIAL_1D_H

#pragma once

#include "DiffusionParameters.h"
#include <cmath>
#include <string>
#include <memory>

namespace special {

/**
 * @class Potential1D
 * @brief Абстрактный базовый класс для одномерного потенциала
 */
class Potential1D {
public:
    virtual ~Potential1D() = default;

    /**
     * @brief Значение потенциала U(x, t)
     */
    virtual double U(double x, double t = 0.0) const = 0;

    /**
     * @brief Сила F(x, t) = -dU/dx
     */
    virtual double F(double x, double t = 0.0) const = 0;

    /**
     * @brief Вторая производная d²U/dx² (для численных схем)
     */
    virtual double U_xx(double x, double t = 0.0) const {
        double h = 1e-6;
        return (U(x + h, t) - 2 * U(x, t) + U(x - h, t)) / (h * h);
    }

    virtual std::string name() const = 0;
};

/**
 * @class FreePotential
 * @brief Свободная диффузия (нет потенциала)
 */
class FreePotential : public Potential1D {
public:
    explicit FreePotential(const DiffusionParameters& params = {})
        : params_(params) {}

    double U(double x, double t = 0.0) const override {
        return 0.0;
    }

    double F(double x, double t = 0.0) const override {
        return 0.0;
    }

    std::string name() const override { return "Free"; }

private:
    DiffusionParameters params_;
};
/**
* @class ConstantForcePotential
* @brief Постоянная сила:  F(x,t) = F = const
*/
class ConstantForcePotential : public Potential1D {
public:
    explicit ConstantForcePotential(double F)
        : F_(F) {}

    double U(double x, double /*t*/) const override {
        return -F_ * x;
    }

    double F(double /*x*/, double /*t*/) const override {
        return F_;
    }

    std::string name() const override {
        return "ConstantForce";
    }

private:
    double F_;
};

/**
 * @class HarmonicPotential
 * @brief Гармонический потенциал U(x) = 0.5*k*x²
 */
class HarmonicPotential : public Potential1D {
public:
    explicit HarmonicPotential(double k = 1.0) : k_(k) {}

    double U(double x, double t = 0.0) const override {
        return 0.5 * k_ * x * x;
    }

    double F(double x, double t = 0.0) const override {
        return -k_ * x;
    }

    std::string name() const override { return "Harmonic"; }

private:
    double k_;
};

/**
 * @class SinPotential
 * @brief Синусоидальный потенциал U(x) = V0 * sin²(πx/L)
 */
class SinPotential : public Potential1D {
public:
    explicit SinPotential(double V0 = 1.0, double L = 1.0)
        : V0_(V0), L_(L) {}

    double U(double x, double t = 0.0) const override {
        double s = std::sin(M_PI * x / L_);
        return V0_ * s * s;
    }

    double F(double x, double t = 0.0) const override {
        // F = -dU/dx = -V0 * 2*sin(πx/L)*cos(πx/L)*(π/L)
        //    = -V0 * sin(2πx/L)*(π/L)
        return -V0_ * std::sin(2.0 * M_PI * x / L_) * (M_PI / L_);
    }

    std::string name() const override { return "Sinusoidal"; }

private:
    double V0_, L_;
};

/**
 * @class RatchetPotential
 * @brief Асимметричный потенциал (рачет)
 * U(x) = V0 * [sin(2πx/L) + 0.25*sin(4πx/L)] / 2
 */
class RatchetPotential : public Potential1D {
public:
    explicit RatchetPotential(double V0 = 1.0, double L = 1.0, double A2 = 0.25)
        : V0_(V0), L_(L), A2_(A2) {}

    double U(double x, double t = 0.0) const override {
        double arg1 = 2.0 * M_PI * x / L_;
        double arg2 = 4.0 * M_PI * x / L_;
        return V0_ * (std::sin(arg1) + A2_ * std::sin(arg2)) / 2.0;
    }

    double F(double x, double t = 0.0) const override {
        // F = -dU/dx
        double arg1 = 2.0 * M_PI * x / L_;
        double arg2 = 4.0 * M_PI * x / L_;
        return -V0_ * (2.0 * M_PI / L_ * std::cos(arg1) +
                       A2_ * 4.0 * M_PI / L_ * std::cos(arg2)) / 2.0;
    }

    std::string name() const override { return "Ratchet"; }

private:
    double V0_, L_, A2_;
};

/**
 * @class TiltingPotential
 * @brief Наклонный потенциал (base + линейная добавка)
 * U_total(x,t) = U_base(x) + F*x + F_osc*sin(ω*t)*x
 */
class TiltingPotential : public Potential1D {
public:
    TiltingPotential(std::unique_ptr<Potential1D> base_potential,
                     double F_static = 0.0,
                     double F_amplitude = 0.0,
                     double omega = 1.0)
        : base_(std::move(base_potential)),
          F_static_(F_static),
          F_amplitude_(F_amplitude),
          omega_(omega) {}

    double U(double x, double t = 0.0) const override {
        double U_base = base_->U(x, t);
        double F_total = F_static_ + F_amplitude_ * std::sin(omega_ * t);
        return U_base - F_total * x;  // Минус, потому что это потенциальная энергия
    }

    double F(double x, double t = 0.0) const override {
        double F_base = base_->F(x, t);
        double F_tilt = F_static_ + F_amplitude_ * std::sin(omega_ * t);
        return F_base + F_tilt;  // Сила = -dU/dx
    }

    std::string name() const override {
        return "Tilting[" + base_->name() + "]";
    }

private:
    std::unique_ptr<Potential1D> base_;
    double F_static_, F_amplitude_, omega_;
};

/**
 * @class FlashingPotential
 * @brief Мигающий потенциал (включается/выключается)
 * U(x,t) = V0*sin²(πx/L) при t_on, 0 при t_off
 */
class FlashingPotential : public Potential1D {
public:
    FlashingPotential(double V0 = 1.0, double L = 1.0,
                      double period = 1.0, double duty_cycle = 0.5)
        : V0_(V0), L_(L), period_(period), duty_cycle_(duty_cycle) {}

    double U(double x, double t = 0.0) const override {
        double phase = std::fmod(t, period_) / period_;
        if (phase < duty_cycle_) {
            double s = std::sin(M_PI * x / L_);
            return V0_ * s * s;
        }
        return 0.0;
    }

    double F(double x, double t = 0.0) const override {
        double phase = std::fmod(t, period_) / period_;
        if (phase < duty_cycle_) {
            return -V0_ * std::sin(2.0 * M_PI * x / L_) * (M_PI / L_);
        }
        return 0.0;
    }

    std::string name() const override { return "Flashing"; }

private:
    double V0_, L_, period_, duty_cycle_;
};

/**
 * @class TravelingWavePotential
 * @brief Движущаяся волна потенциала
 * U(x,t) = V0 * sin²(π(x - vt)/L)
 */
class TravelingWavePotential : public Potential1D {
public:
    TravelingWavePotential(double V0 = 1.0, double L = 1.0, double velocity = 0.1)
        : V0_(V0), L_(L), velocity_(velocity) {}

    double U(double x, double t = 0.0) const override {
        double x_eff = x - velocity_ * t;
        double s = std::sin(M_PI * x_eff / L_);
        return V0_ * s * s;
    }

    double F(double x, double t = 0.0) const override {
        double x_eff = x - velocity_ * t;
        return -V0_ * std::sin(2.0 * M_PI * x_eff / L_) * (M_PI / L_);
    }

    std::string name() const override { return "TravelingWave"; }

private:
    double V0_, L_, velocity_;
};

} // namespace special

#endif // NUMERICAL_METHODS_IN_PHYSICS_POTENTIAL_1D_H