#ifndef NUMERICAL_METHODS_IN_PHYSICS_GAUSSSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_GAUSSSOLVER_H

#pragma once
#include "Base/IIntegralSolver.h"
#include <vector>
#include <string>

// Table values for n = 2,3,4
struct GaussTableEntry {
    std::vector<double> nodes;
    std::vector<double> weights;
};

static GaussTableEntry getGaussTable(size_t n) {
    if (n == 2) { // n = 2
        return {{-0.57735026919, 0.57735026919}, {1.0, 1.0}};
    }
    if (n == 3) {
        return {{-0.77459666924, 0.0, 0.77459666924}, {0.55555555555, 0.88888888889, 0.55555555555}};
    }
    if (n == 4) {
        return {{-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116},
                {0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451}};
    }
    // default: use 2 nodes
    return {{-0.57735026919, 0.57735026919}, {1.0, 1.0}};
}

class GaussSolver : public IIntegralSolver {
public:
    GaussSolver(size_t max_nodes = 4) : max_nodes_(std::max(size_t(2), max_nodes)) {}
    std::string name() const override { return "Gauss-Christoffel"; }

    std::optional<IntegrateResult> integrate(
        std::function<double(double)> func, double a, double b,
        double /*tol*/, size_t /*max_intervals*/ = 0
    ) override {
        std::vector<std::pair<size_t, double>> estim;
        std::vector<double> errors_hist;
        double last_integral = 0.0;

        for (size_t n = 2; n <= max_nodes_; ++n) {
            auto tab = getGaussTable(n);
            double xmid = (a + b) / 2, xhalf = (b - a) / 2;
            double integral = 0.0;
            for (size_t i = 0; i < tab.nodes.size(); ++i)
                integral += tab.weights[i] * func(xmid + xhalf * tab.nodes[i]);
            integral *= xhalf;
            estim.push_back({ n, integral });
            if (n > 2)
                errors_hist.push_back(std::abs(integral - last_integral));
            last_integral = integral;
        }

        double final_result = estim.back().second;
        double final_error = errors_hist.empty() ? 0.0 : errors_hist.back();
        return IntegrateResult{final_result, estim.back().first, final_error, estim, errors_hist};
    }

private:
    size_t max_nodes_;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_GAUSSSOLVER_H