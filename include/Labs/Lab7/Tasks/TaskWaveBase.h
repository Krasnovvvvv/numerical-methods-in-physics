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
#include <functional>

class TaskWaveBase {
public:
    struct PhysParams {
        double c;      // скорость волны, м/с
        double L;      // длина струны, м
        double x0;     // положение удара, м
        double delta;  // полуширина зоны удара, м
        double v0;     // начальная скорость в зоне удара, м/с
    };

    using ExactFunc = std::function<double(double,double,const PhysParams&)>;

    explicit TaskWaveBase(Plotter* plotter = nullptr,
                          ExactFunc exact = nullptr)
        : plotter(plotter)
        , exact(exact)
    {}

    virtual ~TaskWaveBase() = default;

    void run(const PhysParams& phys,
             double tMax_dim, // макс. физ. время, с
             double h_x,      // шаг по x, м
             double tau_t)    // шаг по t, с
    {
        physParams = phys;

        // безразмерные шаги
        double tMax_nd = phys.c * tMax_dim / phys.L;
        double tau_nd  = phys.c * tau_t    / phys.L;
        double h_nd    = h_x / phys.L;

        int N  = static_cast<int>(std::round(1.0 / h_nd));
        int nx = N + 1;
        int M  = static_cast<int>(std::round(tMax_nd / tau_nd));
        int nt = M + 1;

        h_nd   = 1.0 / N;
        tau_nd = tMax_nd / M;

        x.resize(nx);
        t.resize(nt);
        for (int i = 0; i < nx; ++i)
            x[i] = i * h_nd * phys.L;              // x = ξ L
        for (int s = 0; s < nt; ++s)
            t[s] = s * tau_nd * (phys.L / phys.c); // t = τ L/c

        u.resize(nt, nx);
        u.setZero();

        Timer timer;

        solveND(tMax_nd, h_nd, tau_nd, N, M);

        double maxErr = 0.0;
        if (exact)
            maxErr = computeMaxError();

        auto elapsed_us = timer.elapsed();
        std::cout << "--- Wave task --- scheme: " << schemeName()
                  << ", h_x = " << h_nd * phys.L
                  << ", tau_t = " << tau_nd * (phys.L / phys.c)
                  << ", max error = " << maxErr
                  << ", computation time: " << elapsed_us << " microseconds\n";

        if (plotter) {
            plotSnapshots();            // численное решение
            if (exact) {
                plotExact();            // точное решение
                plotAtPointCompare(phys.x0); // сравнение в x0
            }
        }
    }

    void plotAtPoint(double x_fixed) {
        if (!plotter) return;
        int nt_ = static_cast<int>(t.size());
        int nx_ = static_cast<int>(x.size());
        if (nx_ == 0) return;

        int i0 = 0;
        double best = std::abs(x[0] - x_fixed);
        for (int i = 1; i < nx_; ++i) {
            double d = std::abs(x[i] - x_fixed);
            if (d < best) { best = d; i0 = i; }
        }

        std::vector<double> ts(nt_), vals(nt_);
        for (int s = 0; s < nt_; ++s) {
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
    ExactFunc exact;

    Eigen::VectorXd x;      // размерный x
    Eigen::VectorXd t;      // размерный t
    Eigen::MatrixXd u;      // размерное u(x,t)
    PhysParams physParams{};

    virtual void solveND(double tMax_nd,
                         double h_nd,
                         double tau_nd,
                         int N,
                         int M) = 0;

    virtual const char* schemeName() const = 0;

    double computeMaxError() const {
        if (!exact) return 0.0;
        int nx_ = static_cast<int>(x.size());
        int nt_ = static_cast<int>(t.size());
        double maxErr = 0.0;
        for (int s = 0; s < nt_; ++s) {
            for (int i = 0; i < nx_; ++i) {
                double u_ex = exact(x[i], t[s], physParams);
                double diff = std::abs(u(s,i) - u_ex);
                if (diff > maxErr) maxErr = diff;
            }
        }
        return maxErr;
    }

    void plotSnapshots() {
        int nx_ = static_cast<int>(x.size());
        int nt_ = static_cast<int>(t.size());
        if (nx_ == 0 || nt_ == 0) return;

        std::vector<double> rel_times = {0.1, 0.3, 0.5, 1.0};
        double tMax_real = t[nt_ - 1];
        double dt = t[1] - t[0];

        std::vector<double> x_vec(nx_);
        for (int i = 0; i < nx_; ++i) x_vec[i] = x[i];

        std::vector<std::vector<double>> xs, ys;
        std::vector<std::string> labels;

        for (double r : rel_times) {
            double target_t = r * tMax_real;
            int s = static_cast<int>(std::round(target_t / dt));
            if (s >= nt_) s = nt_ - 1;
            double ts = t[s];

            std::vector<double> u_num(nx_);
            for (int i = 0; i < nx_; ++i)
                u_num[i] = u(s, i);

            xs.push_back(x_vec);
            ys.push_back(u_num);
            labels.emplace_back("u_{num}(x,t=" + std::to_string(ts) + ")");
        }

        plotter->plot(xs, ys, labels, "x (m)", "u_{num}(x,t)");
    }

    void plotExact() {
        int nx_ = static_cast<int>(x.size());
        int nt_ = static_cast<int>(t.size());
        if (!exact || nx_ == 0 || nt_ == 0) return;

        std::vector<double> rel_times = {0.1, 0.3, 0.5, 1.0};
        double tMax_real = t[nt_ - 1];
        double dt = t[1] - t[0];

        std::vector<double> x_vec(nx_);
        for (int i = 0; i < nx_; ++i) x_vec[i] = x[i];

        std::vector<std::vector<double>> xs, ys;
        std::vector<std::string> labels;

        for (double r : rel_times) {
            double target_t = r * tMax_real;
            int s = static_cast<int>(std::round(target_t / dt));
            if (s >= nt_) s = nt_ - 1;
            double ts = t[s];

            std::vector<double> u_ex(nx_);
            for (int i = 0; i < nx_; ++i)
                u_ex[i] = exact(x[i], ts, physParams);

            xs.push_back(x_vec);
            ys.push_back(u_ex);
            labels.emplace_back("u_{exact}(x,t=" + std::to_string(ts) + ")");
        }

        plotter->plot(xs, ys, labels, "x (m)", "u_{exact}(x,t)");
    }

    void plotAtPointCompare(double x_fixed) {
        if (!plotter || !exact) return;
        int nt_ = static_cast<int>(t.size());
        int nx_ = static_cast<int>(x.size());
        if (nx_ == 0) return;

        int i0 = 0;
        double best = std::abs(x[0] - x_fixed);
        for (int i = 1; i < nx_; ++i) {
            double d = std::abs(x[i] - x_fixed);
            if (d < best) { best = d; i0 = i; }
        }

        std::vector<double> ts(nt_), u_num(nt_), u_ex(nt_);
        for (int s = 0; s < nt_; ++s) {
            ts[s]   = t[s];
            u_num[s]= u(s, i0);
            u_ex[s] = exact(x_fixed, t[s], physParams);
        }

        std::vector<std::vector<double>> xs = {ts, ts};
        std::vector<std::vector<double>> ys = {u_num, u_ex};
        std::vector<std::string> labels = {"u_{num}(x0,t)","u_{exact}(x0,t)"};

        plotter->plot(xs, ys, labels, "t (s)", "u at x0");
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEBASE_H
