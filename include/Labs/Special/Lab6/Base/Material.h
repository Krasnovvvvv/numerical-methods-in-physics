#ifndef NUMERICAL_METHODS_IN_PHYSICS_MATERIAL_H
#define NUMERICAL_METHODS_IN_PHYSICS_MATERIAL_H
#pragma once

#include "Grid.h"
#include <vector>

namespace fdtd {

    /// Describes a uniform dielectric layer
    struct Layer {
        double x_start;  // physical start position (meters or normalized)
        double x_end;    // physical end position
        double eps_r;    // relative permittivity
        double mu_r  = 1.0;  // relative permeability
        double sigma = 0.0;  // electric conductivity (normalized)
    };

    /// Convenience material presets (non-dispersive refractive indices)
    namespace Materials {
        inline Layer Vacuum(double x_start, double x_end) {
            return {x_start, x_end, 1.0, 1.0, 0.0};
        }
        inline Layer SiO2(double x_start, double x_end) {
            return {x_start, x_end, 1.45 * 1.45, 1.0, 0.0};  // n = 1.45
        }
        inline Layer TiO2(double x_start, double x_end) {
            return {x_start, x_end, 2.28 * 2.28, 1.0, 0.0};  // n = 2.28
        }
        inline Layer Dielectric(double x_start, double x_end, double n) {
            return {x_start, x_end, n * n, 1.0, 0.0};
        }
    }

    /// Apply a set of layers to the grid.
    /// For each layer, sets eps and mu at the corresponding E and H nodes.
    inline void applyLayers(Grid& grid, const std::vector<Layer>& layers) {
        double dx = grid.config().dx;
        for (const auto& layer : layers) {
            // Set eps at E-nodes
            for (int i = 0; i < grid.numE(); ++i) {
                double x = grid.xE(i);
                if (x >= layer.x_start && x < layer.x_end) {
                    grid.eps[i] = layer.eps_r;
                    grid.sigma_e[i] = layer.sigma;
                }
            }
            // Set mu at H-nodes
            for (int i = 0; i < grid.numH(); ++i) {
                double x = grid.xH(i);
                if (x >= layer.x_start && x < layer.x_end) {
                    grid.mu[i] = layer.mu_r;
                }
            }
        }
        grid.recomputeCoefficients();
    }

} // namespace fdtd

#endif //NUMERICAL_METHODS_IN_PHYSICS_MATERIAL_H