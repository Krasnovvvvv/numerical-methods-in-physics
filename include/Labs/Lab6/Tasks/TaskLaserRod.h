#ifndef NUMERICAL_METHODS_IN_PHYSICS_TASKLASERROD_H
#define NUMERICAL_METHODS_IN_PHYSICS_TASKLASERROD_H

#pragma once

#include "Base/ISolver.h"
#include "Helpers/Plotter.h"
#include "Helpers/Timer.h"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

class TaskLaserRod {
public:
    enum class SchemeType {
        ImplicitEuler,   // σ = 1
        CrankNicolson    // σ = 0.5
    };

    TaskLaserRod(ISolver& solver,
                 Plotter* plotter = nullptr,
                 SchemeType scheme = SchemeType::CrankNicolson)
        : solver(solver)
        , plotter(plotter)
        , scheme(scheme)
    {}

    struct PhysParams {
        double rho;  // кг/м^3
        double c;    // Дж/(кг·К)
        double k;    // Вт/(м·К)
        double l;    // м
        double tp;   // с
        double I0;   // Вт/м^2
        double T0;   // К
    };

    void run(const PhysParams& phys,
             double tMax_dim, // макс. физическое время (с)
             double h_x,      // шаг по x (м)
             double tau_t)    // шаг по t (с)
    {
        // 1) обезразмеривание: ξ = x/l, τ = t/tp,
        // Θ = (T - T0)/ΔT*,  ΔT* = I0 tp / k
        double rho = phys.rho;
        double c   = phys.c;
        double k   = phys.k;
        double l   = phys.l;
        double tp  = phys.tp;
        double I0  = phys.I0;

        T0 = phys.T0;
        dT = I0 * tp / k;                          // ΔT*
        double alpha = k * tp / (rho * c * l * l); // коэффициент в ур-нии
        double beta  = I0 / (k * dT);              // = 1 в такой нормировке

        // 2) сетка в ξ, τ
        double h_nd    = h_x      / l;   // шаг по ξ
        double tMax_nd = tMax_dim / tp;  // макс. τ
        double tau_nd  = tau_t    / tp;  // шаг по τ

        const int N  = static_cast<int>(std::round(1.0 / h_nd));     // интервалов по ξ
        const int nx = N + 1;
        const int M  = static_cast<int>(std::round(tMax_nd / tau_nd));
        const int nt = M + 1;

        const double h_eff   = 1.0 * 1.0 / N;      // скорректированный шаг по ξ
        const double tau_eff = tMax_nd / M;        // скорректированный шаг по τ

        const double sigma = (scheme == SchemeType::CrankNicolson) ? 0.5 : 1.0;
        const double mu    = alpha * tau_eff / (h_eff * h_eff);

        // 3) размерные сетки (x,t)
        x.resize(nx);
        t.resize(nt);
        for (int i = 0; i < nx; ++i) x[i] = i * h_eff * l;     // x = ξ l
        for (int s = 0; s < nt; ++s) t[s] = s * tau_eff * tp;  // t = τ tp

        // 4) массив размерной температуры
        u.resize(nt, nx);
        u.setZero();
        for (int i = 0; i < nx; ++i)
            u(0, i) = T0;                 // T(x,0) = T0

        // 5) временной цикл в переменной Θ
        Eigen::MatrixXd theta_layer(nt, nx);
        theta_layer.setZero();            // Θ(ξ,0) = 0

        Timer timer;

        for (int s = 0; s < nt - 1; ++s) {
            Eigen::MatrixXd A(nx, nx);
            Eigen::VectorXd B(nx);
            A.setZero();
            B.setZero();

            // внутренние узлы 1..N-1: θ-схема
            for (int i = 1; i <= N - 1; ++i) {
                int row = i;
                double a_im1 = -sigma * mu;
                double a_i   = 1.0 + 2.0 * sigma * mu;
                double a_ip1 = -sigma * mu;

                if (i - 1 >= 0) A(row, i - 1) += a_im1;
                A(row, i) += a_i;
                if (i + 1 <= N) A(row, i + 1) += a_ip1;

                double uim1  = theta_layer(s, i - 1);
                double ui    = theta_layer(s, i);
                double uip1  = theta_layer(s, i + 1);
                double lap_n = (uip1 - 2.0 * ui + uim1);
                double rhs   = ui + (1.0 - sigma) * mu * lap_n;

                B[row] = rhs;
            }

            // 6) левое граничное условие:
            // ∂Θ/∂ξ|_{0} = -β τ e^{-τ}  =>
            // (Θ1^{s+1} - Θ0^{s+1})/h = -β τ^{s+1} e^{-τ^{s+1}}
            double tau_sp1 = (s + 1) * tau_eff;
            double g_sp1   = beta * tau_sp1 * std::exp(-tau_sp1);


            int rowL = 0;
            A(rowL, 0) = -1.0;
            A(rowL, 1) =  1.0;
            B[rowL]    = -h_eff * g_sp1;

            // 7) правое граничное условие:
            // ∂Θ/∂ξ|_{1} = 0 =>
            // (3Θ_N^{s+1} - 4Θ_{N-1}^{s+1} + Θ_{N-2}^{s+1})/(2h) = 0
            int rowR = N;
            A(rowR, N    ) =  3.0;
            A(rowR, N - 1) = -4.0;
            A(rowR, N - 2) =  1.0;
            B[rowR]       =  0.0;

            auto res = solver.solve(A, B);
            for (int i = 0; i < nx; ++i)
                theta_layer(s + 1, i) = res.solution[i];
        }

        auto elapsed_us = timer.elapsed();
        std::cout << "--- Laser rod solution (dimensional output) ---\n";
        std::cout << "scheme = "
                  << (scheme == SchemeType::CrankNicolson ? "Crank-Nicolson" : "Implicit Euler")
                  << ", h_x = " << h_x
                  << ", tau_t = " << tau_t
                  << ", computation time: " << elapsed_us << " microseconds\n";

        // 8) перевод Θ -> T
        for (int s = 0; s < nt; ++s)
            for (int i = 0; i < nx; ++i)
                u(s, i) = T0 + dT * theta_layer(s, i);

        // 9) размерные графики
        if (plotter) {
            std::vector<double> rel_times = {0.1, 0.3, 0.5, 1.0};
            std::vector<double> x_vec(nx);
            for (int i = 0; i < nx; ++i) x_vec[i] = x[i];

            std::vector<std::vector<double>> xs, ys;
            std::vector<std::string> labels;

            for (double r : rel_times) {
                double target_t = r * tMax_dim;
                int s = static_cast<int>(std::round(target_t / tau_t));
                if (s >= nt) s = nt - 1;
                double ts = t[s];

                std::vector<double> T_num(nx);
                for (int i = 0; i < nx; ++i)
                    T_num[i] = u(s, i);

                xs.push_back(x_vec);
                ys.push_back(T_num);
                labels.emplace_back("T(x,t=" + std::to_string(ts) + ")");
            }

            plotter->plot(xs, ys, labels, "x (m)", "T (K)");
        }
    }

    // T(0,t)
    void plotTemperatureAtLeft() {
        if (!plotter) return;
        int nt_loc = static_cast<int>(t.size());

        std::vector<double> ts(nt_loc), Tleft(nt_loc);
        for (int s = 0; s < nt_loc; ++s) {
            ts[s]    = t[s];
            Tleft[s] = u(s, 0);
        }

        std::vector<std::vector<double>> xs = { ts };
        std::vector<std::vector<double>> ys = { Tleft };
        std::vector<std::string> labels = { "T(0,t)" };

        plotter->plot(xs, ys, labels, "t (s)", "T(0,t) (K)");
    }

    const Eigen::MatrixXd& getU()   const { return u; } // T(x_i,t_s)
    const Eigen::VectorXd& getX()   const { return x; } // x_i (м)
    const Eigen::VectorXd& getTime()const { return t; } // t_s (с)
    double getT0()                  const { return T0; }
    double getDeltaT()              const { return dT; }

private:
    ISolver& solver;
    Plotter* plotter;
    SchemeType scheme;

    // размерные сетки и температура
    Eigen::VectorXd x;  // м
    Eigen::VectorXd t;  // с
    Eigen::MatrixXd u;  // T(x_i,t_s)

    double T0  = 0.0;
    double dT  = 1.0;
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_TASKLASERROD_H
