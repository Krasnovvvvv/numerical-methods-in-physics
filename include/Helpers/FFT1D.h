#ifndef NUMERICAL_METHODS_IN_PHYSICS_FFT1D_H
#define NUMERICAL_METHODS_IN_PHYSICS_FFT1D_H

#pragma once

#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <cmath>

class FFT1D {
public:
    using Complex = std::complex<double>;
    using CVector = Eigen::VectorXcd;
    using DVector = Eigen::VectorXd;

    static CVector forward(const DVector& in) {
        int N = static_cast<int>(in.size());
        CVector a(N);
        for (int i = 0; i < N; ++i) a[i] = Complex(in[i], 0.0);
        fft_recursive(a);
        return a;
    }

    static DVector magnitude(const CVector& spec) {
        int N = static_cast<int>(spec.size());
        DVector mag(N);
        for (int i = 0; i < N; ++i) mag[i] = std::abs(spec[i]);
        return mag;
    }

private:
    static void fft_recursive(CVector& a) {
        int N = static_cast<int>(a.size());
        if (N <= 1) return;

        CVector a_even(N / 2), a_odd(N / 2);
        for (int i = 0; i < N / 2; ++i) {
            a_even[i] = a[2 * i];
            a_odd[i]  = a[2 * i + 1];
        }
        fft_recursive(a_even);
        fft_recursive(a_odd);

        const double PI = std::acos(-1.0);
        for (int k = 0; k < N / 2; ++k) {
            Complex w = std::exp(Complex(0, -2.0 * PI * k / N));
            a[k]         = a_even[k] + w * a_odd[k];
            a[k + N / 2] = a_even[k] - w * a_odd[k];
        }
    }
};

#endif