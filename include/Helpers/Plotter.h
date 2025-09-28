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
    using ExpectedCurveFunc = std::function<std::vector<double>(const std::vector<double> &,const std::vector<double> &)>;

    void plot(const std::vector<double>& x,
              const std::vector<double>& y,
              const std::string& label,
              const std::string& x_label,
              const std::string& y_label,
              bool logarithmic = false,
              std::optional<ExpectedCurveFunc> expected = std::nullopt)
    {
        using namespace matplot;

        if (expected) {
            std::vector<double> y_exp = (*expected)(x,y);

            // Формируем «матрицу» Y (2 строки: экспериментальные и теоретические)
            std::vector<std::vector<double>> ys = {y, y_exp};

            auto lines = matplot::plot(x, ys);

            if(logarithmic)
                lines = matplot::loglog(x,ys);

            lines[0]->marker("o").marker_size(8).color("blue").display_name(label);
            lines[1]->line_style("--").color("red").display_name("Expected");

            double min_y = std::min(*std::min_element(y.begin(), y.end()), *std::min_element(y_exp.begin(), y_exp.end()));
            double max_y = std::max(*std::max_element(y.begin(), y.end()), *std::max_element(y_exp.begin(), y_exp.end()));
            ylim({min_y - 0.1 * std::abs(min_y), max_y + 0.1 * std::abs(max_y)});
        } else {
            auto p = matplot::plot(x, y);
            p->marker_size(12).color("blue").display_name(label);

            double min_y = *std::min_element(y.begin(), y.end());
            double max_y = *std::max_element(y.begin(), y.end());
            ylim({min_y - 0.1 * std::abs(min_y), max_y + 0.1 * std::abs(max_y)});
        }

        xlabel(x_label);
        ylabel(y_label);
        legend()->font_size(10);
        xtickangle(45);
        grid(true);
        show();
    }

};


#endif // NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H
