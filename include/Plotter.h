#ifndef NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H
#define NUMERIAL_METHODS_IN_PHYSICS_PLOTTER_H

#pragma once
#include <matplot/matplot.h>
#include <vector>
#include <string>
#include <functional>
#include <optional>
#include <algorithm>

class Plotter {
public:
    using ExpectedCurveFunc = std::function<std::vector<double>(const std::vector<double>&)>;

    Plotter(std::optional<ExpectedCurveFunc> expected = std::nullopt)
            : expected_curve(expected) {}

    virtual void plot_curve(const std::vector<double>& x, const std::vector<double>& y,
                            const std::string& label) {
        using namespace matplot;

        if (expected_curve) {
            auto y_expected = (*expected_curve)(x);
            auto exp_plot = plot(x, y_expected, "--");
            exp_plot->color("red");
            exp_plot->line_width(2);
            exp_plot->display_name("O(n)");
        }

        auto p = scatter(x, y);
        p->marker_size(8);
        p->color("blue");
        p->display_name(label);

        double min_y = *std::min_element(y.begin(), y.end());
        double max_y = *std::max_element(y.begin(), y.end());
        ylim({min_y - 0.1 * std::abs(min_y), max_y + 0.1 * std::abs(max_y)});

        xlabel("n");
        ylabel("Time, ms");
        legend();
        xtickangle(45);
        show();
    }
private:
    std::optional<ExpectedCurveFunc> expected_curve;
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H
