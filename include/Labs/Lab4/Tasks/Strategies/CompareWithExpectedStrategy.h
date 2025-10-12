#ifndef NUMERICAL_METHODS_IN_PHYSICS_COMPAREWITHEXPECTEDSTRATEGY_H
#define NUMERICAL_METHODS_IN_PHYSICS_COMPAREWITHEXPECTEDSTRATEGY_H

#pragma once
#include <Labs/Lab4/Tasks/IODETaskStrategy.h>

class CompareWithExpectedStrategy : public IODETaskStrategy {
public:
    CompareWithExpectedStrategy(IODESolver* solver,
                               std::function<double(double)> expected_func)
        : solver_(solver), expected_func_(std::move(expected_func)) {}

    void run(const std::function<std::vector<double>(double, const std::vector<double>&)>& rhs,
             const std::vector<double>& y0,
             double t0, double tn,
             double h, double tol,
             unsigned short graphNumber,
             Plotter* plotter,
             std::vector<double> /*steps*/,
             std::vector<double> /*tolerances*/,
             std::vector<IODESolver*> /*solvers*/) override
    {
        if (!solver_) throw std::runtime_error("Solver not set!");
        if (!plotter) throw std::runtime_error("Plotter not set!");
        if (!expected_func_) throw std::runtime_error("Expected function not set!");

        auto num_result = solver_->solve(rhs, y0, t0, tn, h, tol);

        std::vector<double> t = num_result.t;
        std::vector<double> x_num, x_expected, abs_err;
        for (size_t i = 0; i < t.size(); ++i) {
            double xi = num_result.y[i][0];
            double xi_e = expected_func_(t[i]);
            x_num.push_back(xi);
            x_expected.push_back(xi_e);
            abs_err.push_back(std::fabs(xi - xi_e));
        }

        std::vector<std::vector<double>> xs, ys;
        std::vector<std::string> labels;
        if (graphNumber == 1) {
            xs = { t, t };
            ys = { x_num, x_expected };
            labels = { "Numeric", "Expected" };
            plotter->plot(xs, ys, labels, "t", "x");
        } else if (graphNumber == 2) {
            xs = { t };
            ys = { abs_err };
            labels = { "Abs error" };
            plotter->plot(xs, ys, labels, "t", "|x_num - x_expected|");
        }
    }

private:
    IODESolver* solver_;
    std::function<double(double)> expected_func_;
};


#endif //NUMERICAL_METHODS_IN_PHYSICS_COMPAREWITHEXPECTEDSTRATEGY_H