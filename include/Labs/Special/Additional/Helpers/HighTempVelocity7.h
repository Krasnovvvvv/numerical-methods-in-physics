#ifndef NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPVELOCITY7_H
#define NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPVELOCITY7_H
#pragma once
#include "Labs/Special/Additional/Base/RatchetParams.h"
#include <cmath>

inline double highTempVelocity7(const RatchetParams& p) {
    const double pi  = M_PI;
    const double pi2 = pi * pi;

    const double alpha = p.alpha;
    const double V1    = p.V1;
    const double V2    = p.V2;
    const double eps   = p.epsilon;

    const double xi = 1.0 / (4.0 * pi2 * eps);

    const double oneMinusA = 1.0 - alpha;
    const double onePlusA  = 1.0 + alpha;

    const double pref = pi * oneMinusA * oneMinusA * onePlusA * xi;

    const double P51 =
          48.0  * alpha * alpha * std::pow(xi, 4)
        + 220.0 * alpha * alpha * std::pow(xi, 3)
        + 332.0 * alpha * alpha * std::pow(xi, 2)
        + 185.0 * alpha * alpha * xi
        + 34.0  * alpha * alpha
        + 96.0  * alpha * std::pow(xi, 4)
        + 328.0 * alpha * std::pow(xi, 3)
        + 324.0 * alpha * std::pow(xi, 2)
        + 120.0 * alpha * xi
        + 16.0  * alpha
        + 48.0  * std::pow(xi, 4)
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
        + 34.0  * alpha * alpha
        + 256.0 * alpha * std::pow(xi, 6)
        + 1088.0 * alpha * std::pow(xi, 5)
        + 496.0  * alpha * std::pow(xi, 4)
        - 1452.0 * alpha * std::pow(xi, 3)
        - 470.0  * alpha * std::pow(xi, 2)
        + 149.0  * alpha * xi
        - 2.0    * alpha
        + 128.0  * std::pow(xi, 6)
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
               (xi + 4.0) *
               std::pow(1.0 + 4.0 * xi, 4) *
               (4.0 * xi + 9.0));

    const double num71 =
        (9810.0
         + 18.0*alpha*(347.0 + alpha*(718.0 + alpha*(347.0 + 545.0*alpha)))
         + 90081.0*xi
         + 3.0*alpha*(23379.0 + alpha*(40576.0 + 3.0*alpha*(7793.0 + 10009.0*alpha)))*xi
         + (331006.0 + alpha*(332887.0 + alpha*(495410.0 + alpha*(332887.0 + 331006.0*alpha))))*std::pow(xi,2)
         + (627013.0 + alpha*(834040.0 + alpha*(1106866.0 + alpha*(834040.0 + 627013.0*alpha))))*std::pow(xi,3)
         + 2.0*(332489.0 + alpha*(582986.0 + alpha*(731030.0 + alpha*(582986.0 + 332489.0*alpha))))*std::pow(xi,4)
         + 4.0*(101839.0 + alpha*(229296.0 + alpha*(289922.0 + alpha*(229296.0 + 101839.0*alpha))))*std::pow(xi,5)
         + 8.0*(17759.0 + alpha*(49828.0 + alpha*(66074.0 + alpha*(49828.0 + 17759.0*alpha))))*std::pow(xi,6)
         + 128.0*std::pow(1.0 + alpha,2.0)*(203.0 + alpha*(284.0 + 203.0*alpha))*std::pow(xi,7)
         + 1920.0*std::pow(1.0 + alpha,4.0)*std::pow(xi,8));

    const double denom71 =
        512.0 * std::pow(1.0 + xi,3) *
        (4.0 + xi) *
        std::pow(1.0 + 4.0*xi,4) *
        std::pow(9.0 + 4.0*xi,2);

    const double j71 = pref * (num71 / denom71);

    const double num72 =
        (720.0*(149.0 + alpha*(59.0 + alpha*(223.0 + alpha*(59.0 + 149.0*alpha))))
         + 6.0*(1101723.0 + alpha*(407901.0 + alpha*(1226852.0 + 3.0*alpha*(135967.0 + 367241.0*alpha))))*xi
         + (58208581.0 + alpha*(13143379.0 + alpha*(55077542.0 + alpha*(13143379.0 + 58208581.0*alpha))))*std::pow(xi,2)
         + (220124902.0 + alpha*(37296655.0 + alpha*(182630296.0 + alpha*(37296655.0 + 220124902.0*alpha))))*std::pow(xi,3)
         + (440851279.0 + alpha*(94040404.0 + alpha*(320091346.0 + alpha*(94040404.0 + 440851279.0*alpha))))*std::pow(xi,4)
         + 4.0*(126649322.0 + alpha*(44998819.0 + alpha*(80593330.0 + alpha*(44998819.0 + 126649322.0*alpha))))*std::pow(xi,5)
         + 4.0*(87308963.0 + alpha*(50154720.0 + alpha*(53144762.0 + 7.0*alpha*(7164960.0 + 12472709.0*alpha))))*std::pow(xi,6)
         + 16.0*(9184277.0 + alpha*(7936420.0 + alpha*(6529022.0 + alpha*(7936420.0 + 9184277.0*alpha))))*std::pow(xi,7)
         + 512.0*(73001.0 + alpha*(90774.0 + alpha*(75934.0 + alpha*(90774.0 + 73001.0*alpha))))*std::pow(xi,8)
         + 1024.0*(5368.0 + alpha*(9509.0 + alpha*(9418.0 + alpha*(9509.0 + 5368.0*alpha))))*std::pow(xi,9)
         + 1024.0*std::pow(1.0 + alpha,2.0)*(409.0 + alpha*(238.0 + 409.0*alpha))*std::pow(xi,10)
         + 12288.0*std::pow(1.0 + alpha,4.0)*std::pow(xi,11));

    const double denom72 =
        128.0 * std::pow(1.0 + xi,3) *
        std::pow(4.0 + xi,2) *
        std::pow(1.0 + 4.0*xi,5) *
        std::pow(9.0 + 4.0*xi,2) *
        (25.0 + 4.0*xi);

    const double j72 = pref * (num72 / denom72);

    return  j3  * V1 * V1 * V2
          + j51 * std::pow(V1, 4) * V2
          + j52 * V1 * V1 * std::pow(V2, 3)
          + j71 * std::pow(V1, 6) * V2
          + j72 * std::pow(V1, 4) * std::pow(V2, 3);
}

#endif //NUMERICAL_METHODS_IN_PHYSICS_HIGHTEMPVELOCITY7_H