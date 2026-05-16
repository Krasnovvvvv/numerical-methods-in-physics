#ifndef NUMERICAL_METHODS_IN_PHYSICS_TVDLIMITER_H
#define NUMERICAL_METHODS_IN_PHYSICS_TVDLIMITER_H
#pragma once
#include <algorithm>
#include <cmath>

enum class TVDLimiterType {
    VanLeer,
    VanAlbada,
    MinMod,
    Superbee,
    Sweby,
    QUICK,
    UMIST
};

class TVDLimiter {
public:
    static double eval(TVDLimiterType type, double r, double beta = 1.5) {
        if (!std::isfinite(r)) return 0.0;

        switch (type) {
            case TVDLimiterType::VanLeer:
                return (r <= 0.0) ? 0.0 : (2.0 * r) / (1.0 + r);

            case TVDLimiterType::VanAlbada:
                return (r <= 0.0) ? 0.0 : (r * r + r) / (r * r + 1.0);

            case TVDLimiterType::MinMod:
                return std::max(0.0, std::min(1.0, r));

            case TVDLimiterType::Superbee:
                return std::max(0.0, std::max(std::min(2.0 * r, 1.0),
                                              std::min(r, 2.0)));

            case TVDLimiterType::Sweby: {
                beta = std::clamp(beta, 1.0, 2.0);
                return std::max(0.0, std::max(std::min(beta * r, 1.0),
                                              std::min(r, beta)));
            }

            case TVDLimiterType::QUICK:
                return std::max(0.0, std::min({2.0 * r, (3.0 + r) / 4.0, 2.0}));

            case TVDLimiterType::UMIST:
                return std::max(0.0, std::min({
                    2.0 * r,
                    (1.0 + 3.0 * r) / 4.0,
                    (3.0 + r) / 4.0,
                    2.0
                }));
        }

        return 0.0;
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_TVDLIMITER_H