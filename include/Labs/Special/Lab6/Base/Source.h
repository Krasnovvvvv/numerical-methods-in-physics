#pragma once
#include "Grid.h"
#include <cmath>
#include <functional>
#include <stdexcept>

namespace fdtd {

using Waveform = std::function<double(double t)>;

inline Waveform makeCW(double freq, double width = 0.0, double slowness = 3.0) {
    return [=](double t) -> double {
        double envelope = 1.0;
        if (width > 0.0) {
            envelope = 0.5 * (1.0 + std::tanh(t / width - slowness));
        }
        return envelope * std::sin(2.0 * M_PI * freq * t);
    };
}

inline Waveform makeGaussianPulse(double freq_center, double freq_width, double cutoff = 5.0) {
    double w  = 1.0 / freq_width;
    double t0 = cutoff * w;
    return [=](double t) -> double {
        double env = std::exp(-0.5 * (t - t0) * (t - t0) / (w * w));
        return env * std::sin(2.0 * M_PI * freq_center * t);
    };
}

inline Waveform makeGaussianPulseForRange(double lambda_min, double lambda_max,
                                          double c_speed = 1.0, double cutoff = 6.0) {
    double f_max = c_speed / lambda_min;
    double f_min = c_speed / lambda_max;
    double f_center = 0.5 * (f_max + f_min);
    double f_width  = 1.5 * (f_max - f_min);
    return makeGaussianPulse(f_center, f_width, cutoff);
}

enum class SourceMode {
    Soft,       // E += s(t)
    Hard,       // E  = s(t)
    Current     // E += source term equivalent to impressed current density J_y
};

struct Source {
    int         index;
    Waveform    waveform;
    SourceMode  mode = SourceMode::Soft;

    void inject(Grid& grid, int n) const {
        if (!waveform) return;
        if (index < 0 || index >= grid.numE()) {
            throw std::out_of_range("Source index is outside the E-grid");
        }

        double dtp = grid.config().dt_phys();
        double t   = (n + 0.5) * dtp;
        double val = waveform(t);

        switch (mode) {
            case SourceMode::Soft:
                grid.Ey[index] += val;
                break;

            case SourceMode::Hard:
                grid.Ey[index] = val;
                break;

            case SourceMode::Current:
                grid.Ey[index] -= grid.Cb[index] * grid.config().dx * val;
                break;
        }
    }
};

} // namespace fdtd