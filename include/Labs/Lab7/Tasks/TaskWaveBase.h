#ifndef NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEBASE_H
#define NUMERICAL_METHODS_IN_PHYSICS_TASKWAVEBASE_H

#pragma once

#include "Helpers/Plotter.h"
#include "Helpers/Timer.h"

#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

class TaskWaveBase {

public:
    struct PhysParams {
        double c;
        double L;
        double x0;
        double delta;
        double v0;
    };

    using ExactFunc = std::function<double(double, double, const PhysParams&)>;

    enum PlotMode { PLOT_SOLUTION = 1, PLOT_ERROR = 2 };

    explicit TaskWaveBase(Plotter* plotter = nullptr,
                          ExactFunc exact = nullptr,
                          PlotMode mode = PLOT_SOLUTION)
        : plotter(plotter)
        , exact(exact)
        , plotMode(mode)
    {}

    virtual ~TaskWaveBase() = default;

    void run(const PhysParams& phys,
             double tMax,
             double h_x,
             double tau_t)
    {
        physParams = phys;

        // Расчёт в ФИЗИЧЕСКИХ переменных (x в метрах, t в секундах)
        int N  = static_cast<int>(std::round(phys.L / h_x));
        int nx = N + 1;
        int M  = static_cast<int>(std::round(tMax / tau_t));
        int nt = M + 1;

        // Пересчитываем шаги в соответствии с количеством узлов
        h_x   = phys.L / N;
        tau_t = tMax / M;

        x.resize(nx);
        t.resize(nt);

        for (int i = 0; i < nx; ++i)
            x[i] = i * h_x; // координаты в метрах

        for (int s = 0; s < nt; ++s)
            t[s] = s * tau_t; // время в секундах

        u.resize(nt, nx);
        u.setZero();

        Timer timer;

        // Передаём физические шаги
        solveND(tMax, h_x, tau_t, N, M);
        auto elapsed_us = timer.elapsed();

        double maxErr = 0.0;
        if (exact)
            maxErr = computeMaxError();

        std::cout << "--- Wave task --- scheme: " << schemeName()
                  << ", h_x = " << h_x
                  << ", tau_t = " << tau_t
                  << ", max error = " << maxErr
                  << ", time: " << elapsed_us << " us\n";

        if (plotter) {
            if (plotMode == PLOT_SOLUTION)
                plotSnapshots();
            else if (plotMode == PLOT_ERROR && exact)
                plotErrorInTime();
        }
    }

    const Eigen::MatrixXd& getU()   const { return u; }
    const Eigen::VectorXd& getX()   const { return x; }
    const Eigen::VectorXd& getTime() const { return t; }
    const PhysParams&      getPhys() const { return physParams; }

protected:
    Plotter*  plotter;
    ExactFunc exact;
    PlotMode  plotMode;

    Eigen::VectorXd x;
    Eigen::VectorXd t;
    Eigen::MatrixXd u;

    PhysParams physParams{};

    virtual void solveND(double tMax,
                         double h_x,
                         double tau_t,
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
                double diff = std::abs(u(s, i) - u_ex);
                if (diff > maxErr) maxErr = diff;
            }
        }

        return maxErr;
    }

    void plotSnapshots() {
        int nx_ = static_cast<int>(x.size());
        int nt_ = static_cast<int>(t.size());
        if (nx_ == 0 || nt_ == 0) return;

        std::vector<double> rel_times = {0.25, 0.5, 0.75, 1.0};
        double tMax_real = t[nt_ - 1];
        double dt        = (nt_ > 1) ? (t[1] - t[0]) : 1.0;

        std::vector<double> x_vec(nx_);
        for (int i = 0; i < nx_; ++i) x_vec[i] = x[i];

        std::vector<std::vector<double>> xs, ys;
        std::vector<std::string> labels;

        for (double r : rel_times) {
            double target_t = r * tMax_real;
            int s = static_cast<int>(std::round(target_t / dt));
            if (s >= nt_) s = nt_ - 1;
            double ts = t[s];

            // численное решение
            std::vector<double> u_num(nx_);
            for (int i = 0; i < nx_; ++i)
                u_num[i] = u(s, i);

            xs.push_back(x_vec);
            ys.push_back(u_num);
            labels.emplace_back(std::string(schemeName()) +
                                " (t=" + std::to_string(ts).substr(0, 6) + "s)");

            // точное решение, если есть
            if (exact) {
                std::vector<double> u_ex_vec(nx_);
                for (int i = 0; i < nx_; ++i)
                    u_ex_vec[i] = exact(x[i], ts, physParams);

                xs.push_back(x_vec);
                ys.push_back(u_ex_vec);
                labels.emplace_back("Exact (t=" +
                                    std::to_string(ts).substr(0, 6) + "s)");
            }
        }

        if (!xs.empty())
            plotter->plot(xs, ys, labels, "x (m)", "u(x,t)");
    }

    void plotErrorInTime() {
        if (!exact) return;

        int nx_ = static_cast<int>(x.size());
        int nt_ = static_cast<int>(t.size());
        if (nx_ == 0 || nt_ == 0) return;

        std::vector<double> t_vec(nt_);
        std::vector<double> err_vec(nt_);

        for (int s = 0; s < nt_; ++s) {
            t_vec[s] = t[s];
            double maxErr = 0.0;
            for (int i = 0; i < nx_; ++i) {
                double u_ex = exact(x[i], t[s], physParams);
                double diff = std::abs(u(s, i) - u_ex);
                if (diff > maxErr) maxErr = diff;
            }
            err_vec[s] = maxErr;
        }

        std::vector<std::vector<double>> xs(1), ys(1);
        xs[0] = t_vec;
        ys[0] = err_vec;
        std::vector<std::string> labels = { "max |u_{num} - u_{exact}|" };

        plotter->plot(xs, ys, labels, "t (s)", "error");
    }
};

#endif
