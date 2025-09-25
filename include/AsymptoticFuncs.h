#ifndef NUMERICAL_METHODS_IN_PHYSICS_ASYMPTOTICFUNCS_H
#define NUMERICAL_METHODS_IN_PHYSICS_ASYMPTOTICFUNCS_H

#pragma once
#include <vector>

inline std::vector<double> linearAsymptotic(const std::vector<double>& x, const std::vector<double>& y) {
    double sum_x2 = 0, sum_xy = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        sum_x2 += x[i] * x[i];
        sum_xy += x[i] * y[i];
    }
    double alpha = sum_x2 != 0 ? sum_xy / sum_x2 : 0.0;
    std::vector<double> y_expected(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        y_expected[i] = alpha * x[i];
    }
    return y_expected;
}

#endif //NUMERICAL_METHODS_IN_PHYSICS_ASYMPTOTICFUNCS_H
