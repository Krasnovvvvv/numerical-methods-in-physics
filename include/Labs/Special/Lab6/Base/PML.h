#ifndef NUMERICAL_METHODS_IN_PHYSICS_PML_H
#define NUMERICAL_METHODS_IN_PHYSICS_PML_H
#pragma once

#include "Grid.h"
#include <cmath>
#include <algorithm>

namespace fdtd {

/// PML configuration
struct PMLConfig {
    int    thickness = 20;     // number of cells
    int    order     = 3;      // polynomial grading order (0=const, 1=linear, 2=quad, 3=cubic)
    double delta     = 1e-6;   // desired amplitude attenuation
};

/// Apply PML layers on both sides of the grid.
/// Sets sigma_e and sigma_m in the PML regions.
///
/// PML occupies cells [0, thickness) on the left and [numE-thickness, numE) on the right.
///
/// The conductivity profile:  σ(d) = σ_max * (d / L)^m
/// where d = distance from the PML/domain interface, L = PML thickness in physical units.
///
/// Matching condition:  σ* / μ = σ / ε  =>  σ* = σ * μ / ε
/// In free space (ε=μ=1): σ* = σ.
/// More generally, we use the local ε, μ values.
inline void applyPML(Grid& grid, const PMLConfig& cfg) {
    int    N     = cfg.thickness;
    int    m     = cfg.order;
    double dx    = grid.config().dx;
    double L     = N * dx;    // PML physical thickness

    // η₁ = sqrt(μ/ε) at the interface.
    // For simplicity, assume interface is vacuum (η₁ = 1 in normalized units).
    double eta1  = 1.0;

    // σ_max from the desired attenuation
    double sigma_max = -(m + 1.0) * std::log(cfg.delta) / (eta1 * L);

    auto grading = [&](double d) -> double {
        if (d <= 0.0) return 0.0;
        return sigma_max * std::pow(d / L, m);
    };

    int nE = grid.numE();
    int nH = grid.numH();

    // ---- LEFT PML: cells 0..N-1 ----
    // Interface is at cell index N (E-node).
    // Distance for E-node i: d = (N - i) * dx
    for (int i = 0; i < N; ++i) {
        double d = (N - i) * dx;
        double sig = grading(d);
        grid.sigma_e[i] = sig;
    }
    // H-node i corresponds to x_{i+1/2}, distance = (N - i - 0.5)*dx
    for (int i = 0; i < N; ++i) {
        double d = (N - i - 0.5) * dx;
        double sig = grading(d);
        // σ* = σ * μ / ε.  At the PML boundary we use local values (typically 1/1).
        grid.sigma_m[i] = sig * grid.mu[i] / grid.eps[std::min(i, nE - 1)];
    }

    // ---- RIGHT PML: E-nodes [nE-N, nE) ----
    int right_start_E = nE - N;
    for (int i = right_start_E; i < nE; ++i) {
        double d = (i - right_start_E) * dx;
        double sig = grading(d);
        grid.sigma_e[i] = sig;
    }
    // H-nodes: [nH-N, nH)  (approximately)
    int right_start_H = nH - N;
    if (right_start_H < 0) right_start_H = 0;
    for (int i = right_start_H; i < nH; ++i) {
        double d = (i + 0.5 - (nE - N)) * dx;
        if (d < 0) d = 0;
        double sig = grading(d);
        grid.sigma_m[i] = sig * grid.mu[i] / grid.eps[std::min(i, nE - 1)];
    }

    grid.recomputeCoefficients();
}

/// Apply PML only on the right side of the grid.
inline void applyPMLRight(Grid& grid, const PMLConfig& cfg) {
    int    N     = cfg.thickness;
    int    m     = cfg.order;
    double dx    = grid.config().dx;
    double L     = N * dx;
    double eta1  = 1.0;
    double sigma_max = -(m + 1.0) * std::log(cfg.delta) / (eta1 * L);

    auto grading = [&](double d) -> double {
        if (d <= 0.0) return 0.0;
        return sigma_max * std::pow(d / L, m);
    };

    int nE = grid.numE();
    int nH = grid.numH();

    int right_start_E = nE - N;
    for (int i = right_start_E; i < nE; ++i) {
        double d = (i - right_start_E) * dx;
        grid.sigma_e[i] = grading(d);
    }
    int right_start_H = nH - N;
    if (right_start_H < 0) right_start_H = 0;
    for (int i = right_start_H; i < nH; ++i) {
        double d = (i + 0.5 - (nE - N)) * dx;
        if (d < 0) d = 0;
        grid.sigma_m[i] = grading(d) * grid.mu[i] / grid.eps[std::min(i, nE - 1)];
    }
    grid.recomputeCoefficients();
}

/// Apply PML only on the left side of the grid.
inline void applyPMLLeft(Grid& grid, const PMLConfig& cfg) {
    int    N     = cfg.thickness;
    int    m     = cfg.order;
    double dx    = grid.config().dx;
    double L     = N * dx;
    double eta1  = 1.0;
    double sigma_max = -(m + 1.0) * std::log(cfg.delta) / (eta1 * L);

    auto grading = [&](double d) -> double {
        if (d <= 0.0) return 0.0;
        return sigma_max * std::pow(d / L, m);
    };

    int nE = grid.numE();
    int nH = grid.numH();

    for (int i = 0; i < N; ++i) {
        double d = (N - i) * dx;
        grid.sigma_e[i] = grading(d);
    }
    for (int i = 0; i < N; ++i) {
        double d = (N - i - 0.5) * dx;
        grid.sigma_m[i] = grading(d) * grid.mu[i] / grid.eps[std::min(i, nE - 1)];
    }
    grid.recomputeCoefficients();
}

} // namespace fdtd
#endif //NUMERICAL_METHODS_IN_PHYSICS_PML_H