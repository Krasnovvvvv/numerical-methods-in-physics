#ifndef NUMERICAL_METHODS_IN_PHYSICS_ASYMPTOTICFUNCS_H
#define NUMERICAL_METHODS_IN_PHYSICS_ASYMPTOTICFUNCS_H

#pragma once
#include <vector>

inline std::vector<double> linearAsymptotic(const std::vector<double>& x, const std::vector<double>& y) {
    size_t n = x.size();
    double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_xy = 0;
    for (size_t i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_x2 += x[i] * x[i];
        sum_xy += x[i] * y[i];
    }
    double mean_x = sum_x / n;
    double mean_y = sum_y / n;
    double mean_x2 = sum_x2 / n;
    double mean_xy = sum_xy / n;
    double k = (mean_xy - mean_x * mean_y) / (mean_x2 - mean_x * mean_x);
    double b = mean_y - k * mean_x;
    std::vector<double> y_expected(n);
    for (size_t i = 0; i < n; ++i) {
        y_expected[i] = k * x[i] + b;
    }
    return y_expected;
}


#endif //NUMERICAL_METHODS_IN_PHYSICS_ASYMPTOTICFUNCS_H
