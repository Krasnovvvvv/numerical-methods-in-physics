#ifndef NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H
#define NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H

#pragma once
#include <matplot/matplot.h>
#include <vector>
#include <string>
#include <functional>
#include <optional>

class Plotter {
public:
    using ExpectedCurveFunc = std::function<std::vector<double>(const std::vector<double>&)>;

    Plotter(std::optional<ExpectedCurveFunc> expected = std::nullopt)
            : expected_curve(expected) {}

    virtual void plot_curve(const std::vector<double>& x, const std::vector<double>& y,
                            const std::string& label) {
        using namespace matplot;
        auto p = plot(x, y, "-o");
        p->display_name(label);

        if (expected_curve) {
            auto y_expected = (*expected_curve)(x);
            auto exp_plot = plot(x, y_expected, "--");
            exp_plot->display_name("Ожидаемая зависимость");
            exp_plot->color("red");
        }
        xlabel("Размер системы n");
        ylabel("Время, мс");
        legend();
        show();
    }
private:
    std::optional<ExpectedCurveFunc> expected_curve;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H
