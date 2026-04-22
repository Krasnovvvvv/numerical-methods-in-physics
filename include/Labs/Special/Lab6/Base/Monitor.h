#ifndef NUMERICAL_METHODS_IN_PHYSICS_MONITOR_H
#define NUMERICAL_METHODS_IN_PHYSICS_MONITOR_H
#pragma once

#include "Grid.h"
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>

namespace fdtd {

/// DFT monitor at a single spatial point.
/// Accumulates Fourier transform of E_y on-the-fly.
///
/// For 1D normal incidence, we use the amplitude method:
///   R(f) = |E_reflected(f) / E_incident(f)|^2
///   T(f) = |E_transmitted(f) / E_incident(f)|^2
///
/// This is simpler and more robust than the Poynting vector method.
class DFTMonitor {
public:
    DFTMonitor(int index_E, const std::vector<double>& freqs)
        : idx_(index_E), freqs_(freqs)
    {
        int nf = freqs.size();
        E_dft_.assign(nf, {0.0, 0.0});
    }

    /// Accumulate E_y at this monitor point.
    /// Call AFTER E-update and source injection, BEFORE H-update.
    ///
    /// DFT computed as:  Ẽ(f) = Σ_n  E_y(n) · exp(-i·2πf·t_n)
    /// where t_n = (n + 0.5) · dt_phys  is the physical time.
    ///
    /// No dt factor in the sum — we only need amplitude ratios for R/T,
    /// so the DFT values are dimensionless and of order ~1.
    void accumulateAfterE(const Grid& grid, int n) {
        double dtp = grid.config().dt_phys();
        double t   = (n + 0.5) * dtp;       // physical time
        double E_val = grid.Ey[idx_];

        for (int k = 0; k < (int)freqs_.size(); ++k) {
            double omega = 2.0 * M_PI * freqs_[k];
            E_dft_[k] += E_val * std::exp(std::complex<double>(0.0, -omega * t));
        }
    }

    /// Get the complex E-field spectrum
    const std::vector<std::complex<double>>& E_dft() const { return E_dft_; }

    const std::vector<double>& frequencies() const { return freqs_; }
    int index() const { return idx_; }

    void reset() {
        for (auto& v : E_dft_) v = {0.0, 0.0};
    }

private:
    int idx_;
    std::vector<double> freqs_;
    std::vector<std::complex<double>> E_dft_;
};

/// Generate a vector of frequencies for a wavelength range [lambda_min, lambda_max]
inline std::vector<double> freqsForRange(double lambda_min, double lambda_max,
                                          int nf = 500, double c_speed = 1.0) {
    double f_min = c_speed / lambda_max;
    double f_max = c_speed / lambda_min;
    std::vector<double> f(nf);
    for (int i = 0; i < nf; ++i) {
        f[i] = f_min + (f_max - f_min) * i / (nf - 1.0);
    }
    return f;
}

} // namespace fdtd
#endif //NUMERICAL_METHODS_IN_PHYSICS_MONITOR_H