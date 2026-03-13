#ifndef NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPVELOCITY5_H
#define NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPVELOCITY5_H
#pragma once
#include <cmath>

inline double highTempVelocity5(const RatchetParams& p) {
    const double pi = M_PI;
    const double pi2 = pi * pi;

    const double alpha = p.alpha;
    const double V1 = p.V1;
    const double V2 = p.V2;
    const double eps = p.epsilon;

    const double xi = 1.0 / (4.0 * pi2 * eps);

    const double oneMinusA = 1.0 - alpha;
    const double onePlusA  = 1.0 + alpha;
    const double pref = pi * oneMinusA * oneMinusA * onePlusA * xi;

    const double P51 =
          48.0 * alpha * alpha * std::pow(xi, 4)
        + 220.0 * alpha * alpha * std::pow(xi, 3)
        + 332.0 * alpha * alpha * std::pow(xi, 2)
        + 185.0 * alpha * alpha * xi
        + 34.0 * alpha * alpha
        + 96.0 * alpha * std::pow(xi, 4)
        + 328.0 * alpha * std::pow(xi, 3)
        + 324.0 * alpha * std::pow(xi, 2)
        + 120.0 * alpha * xi
        + 16.0 * alpha
        + 48.0 * std::pow(xi, 4)
        + 220.0 * std::pow(xi, 3)
        + 332.0 * std::pow(xi, 2)
        + 185.0 * xi
        + 34.0;

    const double P52 =
          128.0 * alpha * alpha * std::pow(xi, 6)
        + 1568.0 * alpha * alpha * std::pow(xi, 5)
        + 6488.0 * alpha * alpha * std::pow(xi, 4)
        + 9826.0 * alpha * alpha * std::pow(xi, 3)
        + 5455.0 * alpha * alpha * std::pow(xi, 2)
        + 1071.0 * alpha * alpha * xi
        + 34.0 * alpha * alpha
        + 256.0 * alpha * std::pow(xi, 6)
        + 1088.0 * alpha * std::pow(xi, 5)
        + 496.0 * alpha * std::pow(xi, 4)
        - 1452.0 * alpha * std::pow(xi, 3)
        - 470.0 * alpha * std::pow(xi, 2)
        + 149.0 * alpha * xi
        - 2.0 * alpha
        + 128.0 * std::pow(xi, 6)
        + 1568.0 * std::pow(xi, 5)
        + 6488.0 * std::pow(xi, 4)
        + 9826.0 * std::pow(xi, 3)
        + 5455.0 * std::pow(xi, 2)
        + 1071.0 * xi
        + 34.0;

    const double j3 =
        (3.0 * pref * (1.0 + 2.0 * xi)) /
        (4.0 * (1.0 + xi) * std::pow(1.0 + 4.0 * xi, 2));

    const double j51 =
        -(pref * P51) /
        (16.0 * std::pow(1.0 + xi, 2) *
            std::pow(1.0 + 4.0 * xi, 3) *
            (4.0 * xi + 9.0));

    const double j52 =
        -(3.0 * pref * P52) /
        (32.0 * std::pow(1.0 + xi, 2) *
            (xi + 4.0) * std::pow(1.0 + 4.0 * xi, 4) *
            (4.0 * xi + 9.0));

    return j3 * V1 * V1 * V2
         + j51 * V1 * V1 * V1 * V1 * V2
         + j52 * V1 * V1 * V2 * V2 * V2;
}
#endif //NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPVELOCITY5_H