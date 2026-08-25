#ifndef NUMERICAL_METHODS_IN_PHYSICS_LINEARSYSTEM_H
#define NUMERICAL_METHODS_IN_PHYSICS_LINEARSYSTEM_H
#pragma once

#include <vector>
#include <algorithm>
#include "Mesh.h"

struct LinearSystem {
    Mesh mesh;

    std::vector<double> aP, aW, aE, aS, aN, b;

    explicit LinearSystem(const Mesh& m)
        : mesh(m),
          aP(static_cast<size_t>(m.cellCount()), 0.0),
          aW(static_cast<size_t>(m.cellCount()), 0.0),
          aE(static_cast<size_t>(m.cellCount()), 0.0),
          aS(static_cast<size_t>(m.cellCount()), 0.0),
          aN(static_cast<size_t>(m.cellCount()), 0.0),
          b (static_cast<size_t>(m.cellCount()), 0.0) {}

    void clear() {
        std::ranges::fill(aP, 0.0);
        std::ranges::fill(aW, 0.0);
        std::ranges::fill(aE, 0.0);
        std::ranges::fill(aS, 0.0);
        std::ranges::fill(aN, 0.0);
        std::ranges::fill(b,  0.0);
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_LINEARSYSTEM_H