#ifndef NUMERICAL_METHODS_IN_PHYSICS_ADAPTIVESIMPSONSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_ADAPTIVESIMPSONSOLVER_H

#pragma once
#include "Base/IIntegralSolver.h"
#include <cmath>
#include <vector>
#include <functional>
#include <optional>
#include <string>
#include <utility>

class AdaptiveSimpsonSolver : public IIntegralSolver {
public:
    std::string name() const override { return "AdaptiveSimpson"; }

    std::optional<IntegrateResult> integrate(
        std::function<double(double)> func, double a, double b,
        double tol = 1e-8, size_t max_intervals = 1000
    ) override {

        std::vector<std::pair<size_t, double>> estim;
        std::vector<double> errors_hist;
        double last_result = 0;

        for (double curr_tol = 1e-2; curr_tol >= tol; curr_tol /= 2.0) {
            size_t interval_count = 0;
            double last_approx = simpson(func, a, b);
            std::vector<double> tmp_errors;
            double result = adaptiveSimpson(func, a, b,
                curr_tol, max_intervals,
                interval_count,
                last_approx, tmp_errors);
            double error = tmp_errors.empty() ? 0.0 : tmp_errors.back();
            estim.emplace_back(interval_count, result);
            errors_hist.push_back(error);
            last_result = result;
        }

        return IntegrateResult{
            last_result,
            estim.empty() ? 0 : estim.back().first,
            errors_hist.empty() ? 0.0 : errors_hist.back(),
            estim,
            errors_hist
        };
    }

private:
    static double simpson(std::function<double(double)> f, double a, double b) {
        double c = (a + b) / 2.0;
        return (f(a) + 4.0 * f(c) + f(b)) * (b - a) / 6.0;
    }

    static double adaptiveSimpson(
        std::function<double(double)> f, double a, double b, double tol,
        size_t max_intervals, size_t& counter, double S,
        std::vector<double>& errors_hist)
    {
        double c = (a + b) / 2.0;
        double S_left = simpson(f, a, c);
        double S_right = simpson(f, c, b);
        double err = std::abs(S_left + S_right - S) / 15.0;

        if (err < tol || counter + 2 >= max_intervals) {
            counter += 2;
            errors_hist.push_back(err);
            return S_left + S_right;
        } else {
            double left = adaptiveSimpson(f, a, c,
                tol / 2.0, max_intervals, counter, S_left, errors_hist);
            double right = adaptiveSimpson(f, c, b,
                tol / 2.0, max_intervals, counter, S_right, errors_hist);
            return left + right;
        }
    }
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_ADAPTIVESIMPSONSOLVER_H
