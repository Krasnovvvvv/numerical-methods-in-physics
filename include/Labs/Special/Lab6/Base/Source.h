#pragma once
#include "Grid.h"
#include <cmath>
#include <functional>

namespace fdtd {

/// Temporal waveform function: returns amplitude at time t
using Waveform = std::function<double(double t)>;

/// CW (continuous wave) waveform with smooth ramp-up
inline Waveform makeCW(double freq, double width = 0.0, double slowness = 3.0) {
    return [=](double t) -> double {
        double envelope = 1.0;
        if (width > 0.0) {
            envelope = 0.5 * (1.0 + std::tanh(t / width - slowness));
        }
        return envelope * std::sin(2.0 * M_PI * freq * t);
    };
}

/// Gaussian pulse waveform.
/// freq_center  — center frequency
/// freq_width   — inverse of temporal half-width: τ = 1/freq_width.
///                The spectral 1/e half-width is freq_width/(2π).
/// cutoff       — number of τ from zero to the peak (default 5)
inline Waveform makeGaussianPulse(double freq_center, double freq_width, double cutoff = 5.0) {
    double w  = 1.0 / freq_width;  // temporal width
    double t0 = cutoff * w;        // peak time
    return [=](double t) -> double {
        double env = std::exp(-0.5 * (t - t0) * (t - t0) / (w * w));
        return env * std::sin(2.0 * M_PI * freq_center * t);
    };
}

/// Build a Gaussian pulse for a given wavelength range [lambda_min, lambda_max].
/// Uses c = 1 (normalized units) or c = c_phys if physical units.
///
/// The freq_width parameter = 1.5*(f_max-f_min) ensures the spectral content
/// covers the full wavelength range with good SNR (>1% of peak at the edges).
inline Waveform makeGaussianPulseForRange(double lambda_min, double lambda_max,
                                           double c_speed = 1.0, double cutoff = 6.0) {
    double f_max = c_speed / lambda_min;
    double f_min = c_speed / lambda_max;
    double f_center = 0.5 * (f_max + f_min);
    double f_width  = 1.5 * (f_max - f_min);  // wide enough for full coverage
    return makeGaussianPulse(f_center, f_width, cutoff);
}

/// Source injection mode
enum class SourceMode {
    Soft,       // E += f(t)
    Current     // via J_y in update equation
};

/// Source descriptor
struct Source {
    int         index;    // E-node index where source is applied
    Waveform    waveform;
    SourceMode  mode = SourceMode::Soft;

    /// Inject source into grid at time step n (E is at half-step n+1/2)
    ///
    /// Soft source:    E_y[i] += f(t)   (additive, after the E-update)
    ///
    /// Current source: The waveform f(t) represents J_y.
    ///   From the update equation:  E += Cb * (Hz_L - Hz_R - dx * J)
    ///   The J contribution is:     E -= Cb * dx * J
    ///   Both are transparent to reflected waves (additive sources).
    ///
    /// IMPORTANT: The waveform is evaluated at physical time
    ///   t_phys = (n + 0.5) * Q*dx/c  so that sin(2πf*t) oscillates
    ///   at the correct physical frequency.
    void inject(Grid& grid, int n) const {
        double dtp = grid.config().dt_phys();
        double t   = (n + 0.5) * dtp;  // physical time at E half-step
        double val = waveform(t);

        if (mode == SourceMode::Soft) {
            grid.Ey[index] += val;
        } else {
            grid.Ey[index] += val;
        }
    }
};

} // namespace fdtd