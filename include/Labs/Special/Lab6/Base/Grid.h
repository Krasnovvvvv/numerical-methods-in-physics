#ifndef NUMERICAL_METHODS_IN_PHYSICS_GRID_H
#define NUMERICAL_METHODS_IN_PHYSICS_GRID_H
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace fdtd {

/// Parameters for simulation grid
struct GridConfig {
    int    num_cells   = 1000;   // total number of Yee cells
    double dx          = 1.0;    // spatial step (in physical units, e.g. meters)
    double courant     = 0.5;    // Courant number Q = c*dt/dx
    double c_speed     = 1.0;    // speed of light (1.0 for normalized, 3e8 for SI)

    /// Time step used by the solver (normalized units, c=1 effective).
    /// The solver, update coefficients, and PML all use this.
    double dt() const { return courant * dx; }

    /// Physical time step: dt_phys = Q*dx/c.
    /// Used by sources (waveform evaluation) and DFT monitors.
    /// When c_speed=1, dt_phys() == dt().
    double dt_phys() const { return courant * dx / c_speed; }
};

/// 1D Yee grid holding field arrays and material/update coefficients.
/// Convention:
///   Ey[i]    lives at  (x_i, t_{n+1/2})    — integer spatial index
///   Hz[i]    lives at  (x_{i+1/2}, t_n)     — half-integer spatial index
/// Sizes: Ey has num_cells elements,  Hz has num_cells-1 elements.
class Grid {
public:
    explicit Grid(const GridConfig& cfg)
        : cfg_(cfg)
        , n_E_(cfg.num_cells)
        , n_H_(cfg.num_cells - 1)
    {
        Ey  = Eigen::VectorXd::Zero(n_E_);
        Hz  = Eigen::VectorXd::Zero(n_H_);
        eps = Eigen::VectorXd::Ones(n_E_);   // relative permittivity at E nodes
        mu  = Eigen::VectorXd::Ones(n_H_);   // relative permeability at H nodes
        sigma_e = Eigen::VectorXd::Zero(n_E_); // electric conductivity at E nodes
        sigma_m = Eigen::VectorXd::Zero(n_H_); // magnetic conductivity at H nodes

        recomputeCoefficients();
    }

    /// Recompute update coefficients Ca, Cb, Da, Db from current eps, mu, sigma.
    void recomputeCoefficients() {
        double dt = cfg_.dt();
        double dx = cfg_.dx;

        Ca.resize(n_E_);
        Cb.resize(n_E_);
        for (int i = 0; i < n_E_; ++i) {
            double denom = 2.0 * eps[i] + sigma_e[i] * dt;
            Ca[i] = (2.0 * eps[i] - sigma_e[i] * dt) / denom;
            Cb[i] = 2.0 * dt / (dx * denom);
        }

        Da.resize(n_H_);
        Db.resize(n_H_);
        for (int i = 0; i < n_H_; ++i) {
            double denom = 2.0 * mu[i] + sigma_m[i] * dt;
            Da[i] = (2.0 * mu[i] - sigma_m[i] * dt) / denom;
            Db[i] = 2.0 * dt / (dx * denom);
        }
    }

    /// Reset fields to zero (keep material/coefficients)
    void resetFields() {
        Ey.setZero();
        Hz.setZero();
    }

    // Accessors
    int numE() const { return n_E_; }
    int numH() const { return n_H_; }
    const GridConfig& config() const { return cfg_; }

    /// Physical position of E node i
    double xE(int i) const { return i * cfg_.dx; }
    /// Physical position of H node i (at i+1/2)
    double xH(int i) const { return (i + 0.5) * cfg_.dx; }

    // Fields
    Eigen::VectorXd Ey;
    Eigen::VectorXd Hz;

    // Material properties
    Eigen::VectorXd eps;       // ε_r at E-nodes
    Eigen::VectorXd mu;        // μ_r at H-nodes
    Eigen::VectorXd sigma_e;   // σ  at E-nodes
    Eigen::VectorXd sigma_m;   // σ* at H-nodes

    // Update coefficients
    Eigen::VectorXd Ca, Cb;    // for E update
    Eigen::VectorXd Da, Db;    // for H update

private:
    GridConfig cfg_;
    int n_E_, n_H_;
};

} // namespace fdtd

#endif //NUMERICAL_METHODS_IN_PHYSICS_GRID_H