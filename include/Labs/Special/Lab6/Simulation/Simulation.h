#ifndef NUMERICAL_METHODS_IN_PHYSICS_SIMULATION_H
#define NUMERICAL_METHODS_IN_PHYSICS_SIMULATION_H
#pragma once
#include "Labs/Special/Lab6/Base/Grid.h"
#include "Labs/Special/Lab6/Base/Material.h"
#include "Labs/Special/Lab6/Base/Source.h"
#include "Labs/Special/Lab6/Base/PML.h"
#include "Labs/Special/Lab6/Base/Monitor.h"
#include "Labs/Special/Lab6/Solver/FDTDSolver.h"

#include <vector>
#include <functional>
#include <iostream>
#include <string>
#include <utility>
#include <complex>
#include <cmath>
#include <memory>

namespace fdtd {

/// High-level simulation orchestrator.
///
/// Manages grid setup, material/PML/source/monitor configuration, time-stepping,
/// normalization run, and R/T computation using the amplitude method.
///
/// Amplitude method for 1D normal incidence:
///   - Normalization run (no structure): monitor records E_inc(f)
///   - Main run (with structure):
///       reflection monitor records E_total(f) = E_inc(f) + E_ref(f)
///       transmission monitor records E_trans(f)
///   - E_ref(f) = E_total(f) - E_inc(f)
///   - R(f) = |E_ref(f)|² / |E_inc(f)|²
///   - T(f) = |E_trans(f)|² / |E_inc(f)|²
///
class Simulation {
public:
    struct Config {
        int    num_cells      = 1000;
        double dx             = 1.0;
        double courant        = 0.5;
        double c_speed        = 1.0;    // speed of light (1.0=normalized, 3e8=SI)
        int    num_steps      = 10000;
        int    pml_thickness  = 30;
        int    pml_order      = 3;
        double pml_delta      = 1e-6;
    };

    explicit Simulation(const Config& cfg)
        : cfg_(cfg)
    {
        GridConfig gc;
        gc.num_cells = cfg.num_cells;
        gc.dx        = cfg.dx;
        gc.courant   = cfg.courant;
        gc.c_speed   = cfg.c_speed;
        grid_ = std::make_unique<Grid>(gc);
    }

    Grid& grid() { return *grid_; }
    const Grid& grid() const { return *grid_; }

    void addLayer(const Layer& layer) {
        layers_.push_back(layer);
    }

    void setSource(int index, Waveform wf, SourceMode mode = SourceMode::Soft) {
        source_ = Source{index, std::move(wf), mode};
    }

    void addMonitor(int index_E, const std::vector<double>& freqs) {
        monitors_.emplace_back(index_E, freqs);
    }

    /// Apply all configured materials and PML
    void build() {
        if (!layers_.empty()) {
            applyLayers(*grid_, layers_);
        }
        PMLConfig pml_cfg;
        pml_cfg.thickness = cfg_.pml_thickness;
        pml_cfg.order     = cfg_.pml_order;
        pml_cfg.delta     = cfg_.pml_delta;
        applyPML(*grid_, pml_cfg);
    }

    /// Run the simulation
    void run(int snapshot_interval = 0,
             std::function<void(int, const Grid&)> callback = nullptr) {
        for (int n = 0; n < cfg_.num_steps; ++n) {
            FDTDSolver::updateE(*grid_);

            if (source_.waveform) {
                source_.inject(*grid_, n);
            }

            for (auto& mon : monitors_) {
                mon.accumulateAfterE(*grid_, n);
            }

            FDTDSolver::updateH(*grid_);

            if (callback && snapshot_interval > 0 && (n % snapshot_interval == 0)) {
                callback(n, *grid_);
            }
        }
    }

    /// Run normalization (free space + PML, no structure).
    /// Returns monitors with E_inc(f).
    std::vector<DFTMonitor> runNormalization() {
        GridConfig gc;
        gc.num_cells = cfg_.num_cells;
        gc.dx        = cfg_.dx;
        gc.courant   = cfg_.courant;
        gc.c_speed   = cfg_.c_speed;
        Grid norm_grid(gc);

        PMLConfig pml_cfg;
        pml_cfg.thickness = cfg_.pml_thickness;
        pml_cfg.order     = cfg_.pml_order;
        pml_cfg.delta     = cfg_.pml_delta;
        applyPML(norm_grid, pml_cfg);

        std::vector<DFTMonitor> norm_monitors;
        for (const auto& mon : monitors_) {
            norm_monitors.emplace_back(mon.index(), mon.frequencies());
        }

        for (int n = 0; n < cfg_.num_steps; ++n) {
            FDTDSolver::updateE(norm_grid);
            if (source_.waveform) {
                Source norm_src = source_;
                norm_src.inject(norm_grid, n);
            }
            for (auto& mon : norm_monitors) {
                mon.accumulateAfterE(norm_grid, n);
            }
            FDTDSolver::updateH(norm_grid);
        }

        return norm_monitors;
    }

    /// Compute R(f) and T(f) using the amplitude method.
    ///
    /// monitor[0] = reflection point (between source and structure)
    /// monitor[1] = transmission point (after structure)
    ///
    /// norm_monitors = same positions, free-space run
    ///
    /// E_ref(f) = E_total(f) - E_inc(f)     (at reflection point)
    /// E_trans(f) = E_total_trans(f)          (at transmission point)
    /// E_inc(f) = normalization at same point
    ///
    /// R(f) = |E_ref(f)|^2 / |E_inc(f)|^2
    /// T(f) = |E_trans(f)|^2 / |E_inc_trans(f)|^2
    ///
    /// Note: For T in a medium with different refractive index, we'd need
    /// impedance correction. For vacuum→vacuum transmission, amplitude ratio works.

    struct RTResult {
        std::vector<double> freqs;
        std::vector<double> wavelengths;
        std::vector<double> R;
        std::vector<double> T;
    };

    RTResult computeRT(const std::vector<DFTMonitor>& norm_monitors,
                       double c_speed = 1.0) const {
        if (monitors_.size() < 2 || norm_monitors.size() < 2) {
            throw std::runtime_error("Need at least 2 monitors (reflection, transmission)");
        }

        const auto& E_total_ref   = monitors_[0].E_dft();      // total field at refl. point
        const auto& E_total_trans = monitors_[1].E_dft();       // total field at trans. point
        const auto& E_inc_ref     = norm_monitors[0].E_dft();   // incident at refl. point
        const auto& E_inc_trans   = norm_monitors[1].E_dft();   // incident at trans. point

        const auto& freqs = monitors_[0].frequencies();
        int nf = freqs.size();

        // Find max |E_inc|^2 to set a relative threshold.
        // Only compute R/T where the source has meaningful spectral content.
        double max_inc2 = 0;
        for (int k = 0; k < nf; ++k) {
            double e2 = std::norm(E_inc_ref[k]);
            if (e2 > max_inc2) max_inc2 = e2;
        }
        double threshold = max_inc2 * 1e-4;  // 0.01% of peak

        RTResult result;
        result.freqs.resize(nf);
        result.wavelengths.resize(nf);
        result.R.resize(nf);
        result.T.resize(nf);

        for (int k = 0; k < nf; ++k) {
            result.freqs[k] = freqs[k];
            result.wavelengths[k] = c_speed / freqs[k];

            double E_inc_abs2 = std::norm(E_inc_ref[k]);  // |E_inc|^2

            if (E_inc_abs2 > threshold) {
                // Reflected = total - incident
                auto E_ref = E_total_ref[k] - E_inc_ref[k];
                result.R[k] = std::norm(E_ref) / E_inc_abs2;

                // Transmitted (normalized to incident at transmission point)
                double E_inc_trans_abs2 = std::norm(E_inc_trans[k]);
                if (E_inc_trans_abs2 > threshold) {
                    result.T[k] = std::norm(E_total_trans[k]) / E_inc_trans_abs2;
                } else {
                    result.T[k] = std::norm(E_total_trans[k]) / E_inc_abs2;
                }
            } else {
                result.R[k] = 0.0;
                result.T[k] = 0.0;
            }
        }
        return result;
    }

    std::vector<DFTMonitor>& monitors() { return monitors_; }
    const std::vector<DFTMonitor>& monitors() const { return monitors_; }
    const Config& config() const { return cfg_; }

private:
    Config cfg_;
    std::unique_ptr<Grid> grid_;
    std::vector<Layer> layers_;
    Source source_;
    std::vector<DFTMonitor> monitors_;
};

} // namespace fdtd
#endif //NUMERICAL_METHODS_IN_PHYSICS_SIMULATION_H