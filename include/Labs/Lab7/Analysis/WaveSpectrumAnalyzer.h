// WaveSpectrumAnalyzer.h
#ifndef NUMERICAL_METHODS_IN_PHYSICS_WAVESPECTRUMANALYZER_H
#define NUMERICAL_METHODS_IN_PHYSICS_WAVESPECTRUMANALYZER_H

#pragma once

#include "Labs/Lab7/Tasks/TaskWaveBase.h"
#include "Helpers/FFT1D.h"

class WaveSpectrumAnalyzer {
public:
    explicit WaveSpectrumAnalyzer(Plotter* plotter)
        : plotter_(plotter) {}

    // s_index — номер временного слоя, для которого строим спектры
    void analyzeAtTime(const TaskWaveBase& task, int s_index) const {
        const auto& U   = task.getU();
        const auto& x   = task.getX();
        const auto& t   = task.getTime();
        const auto& phys = task.getPhys();

        int nt = static_cast<int>(t.size());
        int nx = static_cast<int>(x.size());
        if (!plotter_ || nx == 0 || nt < 2) return;

        if (s_index < 0) s_index = 0;
        if (s_index >= nt) s_index = nt - 1;

        double h_x = (nx > 1) ? (x[1] - x[0]) : 1.0;
        double tau = (nt > 1) ? (t[1] - t[0]) : 1.0;

        // 1) отклонение u(x, t_s)
        Eigen::VectorXd u_vec(nx);
        for (int i = 0; i < nx; ++i)
            u_vec[i] = U(s_index, i);

        auto U_spec   = FFT1D::forward(u_vec);
        auto U_mag    = FFT1D::magnitude(U_spec);

        // 2) скорость u_t(x, t_s) ~ (u^{s+1}-u^{s-1})/(2τ)
        Eigen::VectorXd ut_vec(nx);
        if (s_index == 0 || s_index == nt - 1) {
            int s1 = (s_index == 0) ? 1 : nt - 2;
            for (int i = 0; i < nx; ++i)
                ut_vec[i] = (U(s1, i) - U(s_index, i)) / tau;
        } else {
            for (int i = 0; i < nx; ++i)
                ut_vec[i] = (U(s_index + 1, i) - U(s_index - 1, i)) / (2.0 * tau);
        }
        auto Ut_spec  = FFT1D::forward(ut_vec);
        auto Ut_mag   = FFT1D::magnitude(Ut_spec);

        // 3) энергия E(x,t_s) = 0.5(ut^2 + c^2 ux^2)
        Eigen::VectorXd ux_vec(nx);
        for (int i = 1; i < nx - 1; ++i)
            ux_vec[i] = (U(s_index, i + 1) - U(s_index, i - 1)) / (2.0 * h_x);
        ux_vec[0]      = (U(s_index, 1) - U(s_index, 0)) / h_x;
        ux_vec[nx - 1] = (U(s_index, nx - 1) - U(s_index, nx - 2)) / h_x;

        Eigen::VectorXd E_vec(nx);
        for (int i = 0; i < nx; ++i)
            E_vec[i] = 0.5 * (ut_vec[i] * ut_vec[i] +
                              phys.c * phys.c * ux_vec[i] * ux_vec[i]);

        auto E_spec  = FFT1D::forward(E_vec);
        auto E_mag   = FFT1D::magnitude(E_spec);

        // ось мод k (0..nx-1)
        std::vector<double> k_vec(nx);
        for (int k = 0; k < nx; ++k) k_vec[k] = k;

        // --- построение трёх спектров ---
        std::vector<std::vector<double>> xs, ys;
        std::vector<std::string> labels;

        xs.push_back(k_vec);
        ys.emplace_back(U_mag.data(), U_mag.data() + nx);
        labels.emplace_back("Spectr u");

        xs.push_back(k_vec);
        ys.emplace_back(Ut_mag.data(), Ut_mag.data() + nx);
        labels.emplace_back("Spectr u_t");

        xs.push_back(k_vec);
        ys.emplace_back(E_mag.data(), E_mag.data() + nx);
        labels.emplace_back("Spectr energy");

        plotter_->plot(xs, ys, labels, "mode k", "amplitude");
    }

private:
    Plotter* plotter_;
};

#endif
