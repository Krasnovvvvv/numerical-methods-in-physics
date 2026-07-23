#ifndef NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H
#define NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H
#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <matplot/matplot.h>
#include "Labs/Special/Lab8/Solver/FlowResult.h"

class Plotter {
public:
    static std::vector<double> yCenters(const FlowResult& r) {
        std::vector<double> y(r.grid.pNy());
        for (int j = 0; j < r.grid.pNy(); ++j) {
            y[j] = (j + 0.5) * r.grid.dy();
        }
        return y;
    }

    static std::vector<double> theoreticalPoiseuille(const FlowResult& r) {
        auto y = yCenters(r);
        std::vector<double> u(y.size());
        for (size_t k = 0; k < y.size(); ++k) {
            const double eta = y[k] / r.cfg.Ly;
            u[k] = 6.0 * r.cfg.uMean * eta * (1.0 - eta);
        }
        return u;
    }

    static std::pair<std::vector<double>, std::vector<double>> profileAtX(const FlowResult& r, double x) {
        int i = std::clamp(static_cast<int>(x / r.grid.dx()), 0, r.grid.pNx() - 1);
        std::vector<double> y(r.grid.pNy());
        std::vector<double> u(r.grid.pNy());
        for (int j = 0; j < r.grid.pNy(); ++j) {
            y[j] = (j + 0.5) * r.grid.dy();
            u[j] = 0.5 * (r.u(i, j) + r.u(i + 1, j));
        }
        return {y, u};
    }

    static std::pair<std::vector<double>, std::vector<double>> centerlinePressure(const FlowResult& r) {
        const int j = r.grid.pNy() / 2;
        std::vector<double> x(r.grid.pNx());
        std::vector<double> p(r.grid.pNx());
        for (int i = 0; i < r.grid.pNx(); ++i) {
            x[i] = (i + 0.5) * r.grid.dx();
            p[i] = r.p(i, j);
        }
        return {x, p};
    }

    static std::pair<std::vector<double>, std::vector<double>> centerlineU(const FlowResult& r) {
        const int j = r.grid.pNy() / 2;
        std::vector<double> x(r.grid.pNx());
        std::vector<double> u(r.grid.pNx());
        for (int i = 0; i < r.grid.pNx(); ++i) {
            x[i] = (i + 0.5) * r.grid.dx();
            u[i] = 0.5 * (r.u(i, j) + r.u(i + 1, j));
        }
        return {x, u};
    }

    static std::vector<std::vector<double>> pressureMatrix(const FlowResult& r) {
        std::vector<std::vector<double>> m(r.grid.pNy(), std::vector<double>(r.grid.pNx()));
        for (int j = 0; j < r.grid.pNy(); ++j) {
            for (int i = 0; i < r.grid.pNx(); ++i) {
                m[r.grid.pNy() - 1 - j][i] = r.p(i, j);
            }
        }
        return m;
    }

    static std::vector<std::vector<double>> uMatrix(const FlowResult& r) {
        std::vector<std::vector<double>> m(r.grid.pNy(), std::vector<double>(r.grid.pNx()));
        for (int j = 0; j < r.grid.pNy(); ++j) {
            for (int i = 0; i < r.grid.pNx(); ++i) {
                m[r.grid.pNy() - 1 - j][i] = 0.5 * (r.u(i, j) + r.u(i + 1, j));
            }
        }
        return m;
    }

    static std::vector<std::vector<double>> vMatrix(const FlowResult& r) {
        std::vector<std::vector<double>> m(r.grid.pNy(), std::vector<double>(r.grid.pNx()));
        for (int j = 0; j < r.grid.pNy(); ++j) {
            for (int i = 0; i < r.grid.pNx(); ++i) {
                m[r.grid.pNy() - 1 - j][i] = 0.5 * (r.v(i, j) + r.v(i, j + 1));
            }
        }
        return m;
    }

    static void plotFieldOverview(const FlowResult& r, const std::string& prefix) {
        using namespace matplot;

        figure(true);
        hold(on);
        imagesc(uMatrix(r));
        colorbar();
        title(prefix + " : u(x,y)");
        xlabel("i");
        ylabel("j");
        show();

        figure(true);
        imagesc(vMatrix(r));
        colorbar();
        title(prefix + " : v(x,y)");
        xlabel("i");
        ylabel("j");
        show();

        figure(true);
        imagesc(pressureMatrix(r));
        colorbar();
        title(prefix + " : p(x,y)");
        xlabel("i");
        ylabel("j");

        show();
        hold(off);
    }

    static void plotInletComparison(const FlowResult& uniformCase,
                                    const FlowResult& profileCase) {
        using namespace matplot;
        auto [y1, u1] = profileAtX(uniformCase, 0.85 * uniformCase.grid.Lx);
        auto [y2, u2] = profileAtX(profileCase, 0.85 * profileCase.grid.Lx);
        auto uth = theoreticalPoiseuille(profileCase);

        figure(true);
        auto l1 = plot(u1, y1, "-o");
        l1 -> line_width(2);
        hold(on);
        auto l2 = plot(u2, y2, "-s");
        l2 -> line_width(2);
        auto l3 = plot(uth, y2, "--");
        l3 -> line_width(2);
        hold(off);
        xlabel("u");
        ylabel("y");
        title("Сравнение профиля скорости: uniform inlet vs parabolic inlet");
        legend({"uniform inlet", "parabolic inlet", "theoretical Poiseuille"});
        show();
    }

    static void plotOutletStudy(const std::vector<double>& lengths,
                                const std::vector<FlowResult>& results) {
        using namespace matplot;
        figure(true);
        hold(on);
        for (size_t k = 0; k < results.size(); ++k) {
            auto [y, u] = profileAtX(results[k], 0.85 * results[k].grid.Lx);
            plot(u, y);
        }
        const auto uth = theoreticalPoiseuille(results.back());
        const auto y = yCenters(results.back());
        auto l1 = plot(uth, y, "k--");
        l1 -> line_width(2);
        hold(off);
        xlabel("u");
        ylabel("y");
        title("Влияние положения outlet на профиль скорости");

        std::vector<std::string> labels;
        for (double l : lengths) {
            labels.push_back("L=" + std::to_string(l));
        }
        labels.push_back("Poiseuille");
        legend(labels);
        show();

        figure(true);
        hold(on);
        for (size_t k = 0; k < results.size(); ++k) {
            auto [x, p] = centerlinePressure(results[k]);
            auto l = plot(x, p);
            l -> line_width(2);
        }
        hold(off);
        xlabel("x");
        ylabel("p centerline");
        title("Сравнение давления при разных outlet положениях");
        legend({"L=" + std::to_string(lengths[0]),
                "L=" + std::to_string(lengths[1]),
                "L=" + std::to_string(lengths[2])});
        show();
    }

    static void plotReStudy(const std::vector<double>& re,
                            const std::vector<double>& yPlus) {
        using namespace matplot;
        figure(true);
        auto l1 = plot(re, yPlus, "-o");
        l1 -> line_width(2);
        hold(on);
        auto l2 = plot(std::vector<double>{re.front(), re.back()}, std::vector<double>{11.63, 11.63}, "r--");
        hold(off);
        l2 -> line_width(2);
        xlabel("Re");
        ylabel("max y+");
        title("Рост y+ при увеличении Re");
        legend({"computed max y+", "y+=11.63"});
        show();
    }
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H