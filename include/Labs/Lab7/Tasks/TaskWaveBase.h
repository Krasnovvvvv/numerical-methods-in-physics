#ifndef NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEBASE_H
#define NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEBASE_H

#pragma once

#include "Helpers/Plotter.h"
#include "Helpers/Timer.h"

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

class TaskWaveBase {
public:
    struct PhysParams {
        double c;      // скорость волны, м/с
        double L;      // длина струны, м
        double x0;     // положение удара, м
        double delta;  // полуширина зоны удара, м
        double v0;     // скорость, сообщённая в зоне удара, м/с
    };

    explicit TaskWaveBase(Plotter* plotter = nullptr)
        : plotter(plotter)
    {}

    virtual ~TaskWaveBase() = default;

    // единый интерфейс запуска для всех схем
    void run(const PhysParams& phys,
             double tMax_dim, // макс. физ. время (с)
             double h_x,      // шаг по x (м)
             double tau_t)    // шаг по t (с)
    {
        physParams = phys;

        // безразмерные шаги
        double tMax_nd = phys.c * tMax_dim / phys.L;
        double tau_nd  = phys.c * tau_t    / phys.L;
        double h_nd    = h_x / phys.L;

        // грубые N, M
        int N  = static_cast<int>(std::round(1.0 / h_nd));
        int nx = N + 1;
        int M  = static_cast<int>(std::round(tMax_nd / tau_nd));
        int nt = M + 1;

        // фактические безразмерные шаги
        h_nd   = 1.0 / N;
        tau_nd = tMax_nd / M;

        // сохраняем размерную сетку
        x.resize(nx);
        t.resize(nt);
        for (int i = 0; i < nx; ++i)
            x[i] = i * h_nd * phys.L;            // x = ξ L
        for (int s = 0; s < nt; ++s)
            t[s] = s * tau_nd * (phys.L / phys.c); // t = τ L/c

        u.resize(nt, nx);
        u.setZero();

        Timer timer;

        // реальное численное решение делегируем наследнику
        solveND(tMax_nd, h_nd, tau_nd, N, M);

        auto elapsed_us = timer.elapsed();
        std::cout << "--- Wave task --- scheme: " << schemeName()
                  << ", h_x = " << h_nd * phys.L
                  << ", tau_t = " << tau_nd * (phys.L / phys.c)
                  << ", computation time: " << elapsed_us << " microseconds\n";

        if (plotter)
            plotSnapshots();
    }

    // u(x_fixed, t)
    void plotAtPoint(double x_fixed) {
        if (!plotter) return;
        int nt = static_cast<int>(t.size());
        int nx = static_cast<int>(x.size());
        if (nx == 0) return;

        int i0 = 0;
        double best = std::abs(x[0] - x_fixed);
        for (int i = 1; i < nx; ++i) {
            double d = std::abs(x[i] - x_fixed);
            if (d < best) { best = d; i0 = i; }
        }

        std::vector<double> ts(nt), vals(nt);
        for (int s = 0; s < nt; ++s) {
            ts[s]   = t[s];
            vals[s] = u(s, i0);
        }

        std::vector<std::vector<double>> xs = {ts};
        std::vector<std::vector<double>> ys = {vals};
        std::vector<std::string> labels = {"u(x=" + std::to_string(x[i0]) + ",t)"};
        plotter->plot(xs, ys, labels, "t (s)", "u");
    }

    const Eigen::MatrixXd& getU()    const { return u; }
    const Eigen::VectorXd& getX()    const { return x; }
    const Eigen::VectorXd& getTime() const { return t; }
    const PhysParams&      getPhys() const { return physParams; }

protected:
    Plotter* plotter;
    Eigen::VectorXd x;      // размерный x
    Eigen::VectorXd t;      // размерный t
    Eigen::MatrixXd u;      // размерное u(x,t)
    PhysParams physParams{};

    // реализация конкретной схемы в безразмерных переменных:
    // θ_{ττ} = θ_{ξξ}, ξ∈[0,1], τ∈[0,tMax_nd]
    virtual void solveND(double tMax_nd,
                         double h_nd,
                         double tau_nd,
                         int N,
                         int M) = 0;

    virtual const char* schemeName() const = 0;

    // отрисовка u(x,t) в несколько моментов времени
    void plotSnapshots() {
        int nx = static_cast<int>(x.size());
        int nt = static_cast<int>(t.size());
        if (nx == 0 || nt == 0) return;

        std::vector<double> rel_times = {0.1, 0.3, 0.5, 1.0};
        double tMax_real = t[nt - 1];
        double dt = t[1] - t[0];

        std::vector<double> x_vec(nx);
        for (int i = 0; i < nx; ++i) x_vec[i] = x[i];

        std::vector<std::vector<double>> xs, ys;
        std::vector<std::string> labels;

        for (double r : rel_times) {
            double target_t = r * tMax_real;
            int s = static_cast<int>(std::round(target_t / dt));
            if (s >= nt) s = nt - 1;
            double ts = t[s];

            std::vector<double> u_num(nx);
            for (int i = 0; i < nx; ++i)
                u_num[i] = u(s, i);

            xs.push_back(x_vec);
            ys.push_back(u_num);
            labels.emplace_back("u(x,t=" + std::to_string(ts) + ")");
        }

        plotter->plot(xs, ys, labels, "x (m)", "u(x,t)");
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEBASE_H
