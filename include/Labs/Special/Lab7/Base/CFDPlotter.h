#ifndef NUMERICAL_METHODS_IN_PHYSICS_CFDPLOTTER_H
#define NUMERICAL_METHODS_IN_PHYSICS_CFDPLOTTER_H
#pragma once

#include <matplot/matplot.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <utility>

#include "Mesh.h"
#include "ScalarField.h"

class CFDPlotter {
public:
    static std::vector<double> xCenters(const Mesh& mesh) {
        std::vector<double> x(static_cast<size_t>(mesh.Nx));
        for (int i = 0; i < mesh.Nx; ++i) {
            x[static_cast<size_t>(i)] = mesh.xc(i);
        }
        return x;
    }

    static std::vector<double> yCenters(const Mesh& mesh) {
        std::vector<double> y(static_cast<size_t>(mesh.Ny));
        for (int j = 0; j < mesh.Ny; ++j) {
            y[static_cast<size_t>(j)] = mesh.yc(j);
        }
        return y;
    }

    static std::vector<std::vector<double>> fieldMatrix(const ScalarField& phi) {
        std::vector<std::vector<double>> z(
            static_cast<size_t>(phi.mesh.Ny),
            std::vector<double>(static_cast<size_t>(phi.mesh.Nx), 0.0)
        );

        for (int j = 0; j < phi.mesh.Ny; ++j) {
            for (int i = 0; i < phi.mesh.Nx; ++i) {
                z[static_cast<size_t>(j)][static_cast<size_t>(i)] = phi(i, j);
            }
        }

        return z;
    }

    static std::pair<std::vector<double>, std::vector<double>>
    extractProfileX(const ScalarField& phi, double y0) {
        const int j = closestJ(phi.mesh, y0);

        std::vector<double> x = xCenters(phi.mesh);
        std::vector<double> c(static_cast<size_t>(phi.mesh.Nx), 0.0);

        for (int i = 0; i < phi.mesh.Nx; ++i) {
            c[static_cast<size_t>(i)] = phi(i, j);
        }

        return {x, c};
    }

    static std::pair<std::vector<double>, std::vector<double>>
    extractProfileY(const ScalarField& phi, double x0) {
        const int i = closestI(phi.mesh, x0);

        std::vector<double> y = yCenters(phi.mesh);
        std::vector<double> c(static_cast<size_t>(phi.mesh.Ny), 0.0);

        for (int j = 0; j < phi.mesh.Ny; ++j) {
            c[static_cast<size_t>(j)] = phi(i, j);
        }

        return {y, c};
    }

    static double sampleBilinear(const ScalarField& phi, double x, double y) {
        const Mesh& m = phi.mesh;

        const double xMin = m.xc(0);
        const double xMax = m.xc(m.Nx - 1);
        const double yMin = m.yc(0);
        const double yMax = m.yc(m.Ny - 1);

        x = std::clamp(x, xMin, xMax);
        y = std::clamp(y, yMin, yMax);

        if (m.Nx == 1 && m.Ny == 1) {
            return phi(0, 0);
        }

        if (m.Nx == 1) {
            int j0 = static_cast<int>(std::floor((y - yMin) / m.dy));
            j0 = std::clamp(j0, 0, m.Ny - 2);
            const int j1 = j0 + 1;
            const double y0 = m.yc(j0);
            const double y1 = m.yc(j1);
            const double ty = (std::abs(y1 - y0) < 1e-14) ? 0.0 : (y - y0) / (y1 - y0);
            return (1.0 - ty) * phi(0, j0) + ty * phi(0, j1);
        }

        if (m.Ny == 1) {
            int i0 = static_cast<int>(std::floor((x - xMin) / m.dx));
            i0 = std::clamp(i0, 0, m.Nx - 2);
            const int i1 = i0 + 1;
            const double x0 = m.xc(i0);
            const double x1 = m.xc(i1);
            const double tx = (std::abs(x1 - x0) < 1e-14) ? 0.0 : (x - x0) / (x1 - x0);
            return (1.0 - tx) * phi(i0, 0) + tx * phi(i1, 0);
        }

        int i0 = static_cast<int>(std::floor((x - xMin) / m.dx));
        int j0 = static_cast<int>(std::floor((y - yMin) / m.dy));

        i0 = std::clamp(i0, 0, m.Nx - 2);
        j0 = std::clamp(j0, 0, m.Ny - 2);

        const int i1 = i0 + 1;
        const int j1 = j0 + 1;

        const double x0 = m.xc(i0);
        const double x1 = m.xc(i1);
        const double y0 = m.yc(j0);
        const double y1 = m.yc(j1);

        const double tx = (std::abs(x1 - x0) < 1e-14) ? 0.0 : (x - x0) / (x1 - x0);
        const double ty = (std::abs(y1 - y0) < 1e-14) ? 0.0 : (y - y0) / (y1 - y0);

        const double c00 = phi(i0, j0);
        const double c10 = phi(i1, j0);
        const double c01 = phi(i0, j1);
        const double c11 = phi(i1, j1);

        const double c0 = (1.0 - tx) * c00 + tx * c10;
        const double c1 = (1.0 - tx) * c01 + tx * c11;

        return (1.0 - ty) * c0 + ty * c1;
    }

    static std::pair<std::vector<double>, std::vector<double>>
    sampleSegment(const ScalarField& phi,
                  double x0, double y0,
                  double x1, double y1,
                  size_t points = 401) {
        if (points < 2) {
            throw std::runtime_error("sampleSegment: points must be >= 2.");
        }

        std::vector<double> s(points, 0.0);
        std::vector<double> c(points, 0.0);

        for (size_t k = 0; k < points; ++k) {
            const double t = static_cast<double>(k) / static_cast<double>(points - 1);
            const double x = (1.0 - t) * x0 + t * x1;
            const double y = (1.0 - t) * y0 + t * y1;

            s[k] = t;
            c[k] = sampleBilinear(phi, x, y);
        }

        return {s, c};
    }

    static void plotHeatmap(const ScalarField& phi,
                        const std::string& titleText = "c(x,y)") {
        using namespace matplot;

        std::vector<std::vector<double>> X(
            static_cast<size_t>(phi.mesh.Ny),
            std::vector<double>(static_cast<size_t>(phi.mesh.Nx), 0.0)
        );

        std::vector<std::vector<double>> Y(
            static_cast<size_t>(phi.mesh.Ny),
            std::vector<double>(static_cast<size_t>(phi.mesh.Nx), 0.0)
        );

        std::vector<std::vector<double>> Z(
            static_cast<size_t>(phi.mesh.Ny),
            std::vector<double>(static_cast<size_t>(phi.mesh.Nx), 0.0)
        );

        for (int j = 0; j < phi.mesh.Ny; ++j) {
            for (int i = 0; i < phi.mesh.Nx; ++i) {
                const size_t row = static_cast<size_t>(j);
                const size_t col = static_cast<size_t>(i);

                X[row][col] = phi.mesh.xc(i);
                Y[row][col] = phi.mesh.yc(j);
                Z[row][col] = phi(i, j);
            }
        }

        figure();

        auto c = contourf(X, Y, Z);

        colorbar();
        colormap(palette::viridis());

        xlabel("x");
        ylabel("y");
        title(titleText);

        xlim({0.0, phi.mesh.Lx});
        ylim({0.0, phi.mesh.Ly});

        grid(off);
        box(on);

        show();
    }

    static void plotProfileX(const ScalarField& phi,
                             double y0,
                             const std::string& label = "c(x, y0)",
                             const std::string& titleText = "Profile along x") {
        using namespace matplot;

        auto [x, c] = extractProfileX(phi, y0);

        figure();
        auto p = plot(x, c);
        p->line_width(2);
        p->display_name(label);
        xlabel("x");
        ylabel("c");
        title(titleText);
        legend();
        grid(on);
        show();
    }

    static void plotProfileY(const ScalarField& phi,
                             double x0,
                             const std::string& label = "c(x0, y)",
                             const std::string& titleText = "Profile along y") {
        using namespace matplot;

        auto [y, c] = extractProfileY(phi, x0);

        figure();
        auto p = plot(y, c);
        p->line_width(2);
        p->display_name(label);
        xlabel("y");
        ylabel("c");
        title(titleText);
        legend();
        grid(on);
        show();
    }

    static void plotCurves(const std::vector<double>& x,
                           const std::vector<std::vector<double>>& ys,
                           const std::vector<std::string>& labels,
                           const std::string& titleText,
                           const std::string& xLabel,
                           const std::string& yLabel) {
        using namespace matplot;

        if (ys.empty()) {
            throw std::runtime_error("plotCurves: empty curve set.");
        }
        if (labels.size() != ys.size()) {
            throw std::runtime_error("plotCurves: labels size mismatch.");
        }

        figure();
        hold(on);

        for (size_t k = 0; k < ys.size(); ++k) {
            auto p = plot(x, ys[k]);
            p->line_width(2);
            p->display_name(labels[k]);
        }

        xlabel(xLabel);
        ylabel(yLabel);
        title(titleText);
        legend();
        grid(on);
        hold(off);
        show();
    }

    static void plotHistory(const std::vector<double>& history,
                            const std::string& titleText = "Outer iteration history") {
        using namespace matplot;

        if (history.empty()) {
            return;
        }

        std::vector<double> it(history.size(), 0.0);
        for (size_t k = 0; k < history.size(); ++k) {
            it[k] = static_cast<double>(k + 1);
        }

        figure();
        auto p = semilogy(it, history);
        p->line_width(2);
        p->marker("o");
        xlabel("Outer iteration");
        ylabel("max|c^(k+1)-c^k|");
        title(titleText);
        grid(on);
        show();
    }

private:
    static int closestI(const Mesh& mesh, double x0) {
        int best = 0;
        double bestDist = std::abs(mesh.xc(0) - x0);

        for (int i = 1; i < mesh.Nx; ++i) {
            const double d = std::abs(mesh.xc(i) - x0);
            if (d < bestDist) {
                bestDist = d;
                best = i;
            }
        }

        return best;
    }

    static int closestJ(const Mesh& mesh, double y0) {
        int best = 0;
        double bestDist = std::abs(mesh.yc(0) - y0);

        for (int j = 1; j < mesh.Ny; ++j) {
            const double d = std::abs(mesh.yc(j) - y0);
            if (d < bestDist) {
                bestDist = d;
                best = j;
            }
        }

        return best;
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_CFDPLOTTER_H