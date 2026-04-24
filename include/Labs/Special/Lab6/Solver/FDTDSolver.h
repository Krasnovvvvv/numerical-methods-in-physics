#ifndef NUMERICAL_METHODS_IN_PHYSICS_FDTDSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_FDTDSOLVER_H
#pragma once

#include "Labs/Special/Lab6/Base/Grid.h"

namespace fdtd {

    /// Perform one leap-frog step of the 1D FDTD algorithm.
    ///
    /// Call order within one time step n:
    ///   1. updateE(grid)           — advances E from n-1/2 to n+1/2
    ///   2. source.inject(grid, n)  — inject source into E at n+1/2
    ///   3. monitors.accumulate()   — record E at n+1/2, H at n
    ///   4. updateH(grid)           — advances H from n to n+1
    ///
    class FDTDSolver {
    public:
        /// Update E_y: n-1/2 → n+1/2
        /// E_y^{n+1/2}_i = Ca_i * E_y^{n-1/2}_i + Cb_i * (H_z^n_{i-1/2} - H_z^n_{i+1/2})
        static void updateE(Grid& grid) {
            int nE = grid.numE();
            for (int i = 1; i < nE - 1; ++i) {
                grid.Ey[i] = grid.Ca[i] * grid.Ey[i]
                            + grid.Cb[i] * (grid.Hz[i - 1] - grid.Hz[i]);
            }
            // Boundary E-nodes: if no PML wraps them, they remain zero (PEC)
            // or are handled by PML loss terms.
            // We still update them if they have valid neighbors:
            // i=0: only has Hz[0] on the right
            grid.Ey[0] = grid.Ca[0] * grid.Ey[0]
                        + grid.Cb[0] * (0.0 - grid.Hz[0]);  // Hz[-1/2] = 0 (PEC)
            // i=nE-1: only has Hz[nE-2] on the left
            grid.Ey[nE - 1] = grid.Ca[nE - 1] * grid.Ey[nE - 1]
                             + grid.Cb[nE - 1] * (grid.Hz[nE - 2] - 0.0);  // PEC right
        }

        /// Update H_z: n → n+1
        /// H_z^{n+1}_{i+1/2} = Da_{i+1/2} * H_z^n_{i+1/2} + Db_{i+1/2} * (E_y^{n+1/2}_i - E_y^{n+1/2}_{i+1})
        static void updateH(Grid& grid) {
            int nH = grid.numH();
            for (int i = 0; i < nH; ++i) {
                grid.Hz[i] = grid.Da[i] * grid.Hz[i]
                            + grid.Db[i] * (grid.Ey[i] - grid.Ey[i + 1]);
            }
        }
    };

} // namespace fdtd
#endif //NUMERICAL_METHODS_IN_PHYSICS_FDTDSOLVER_H