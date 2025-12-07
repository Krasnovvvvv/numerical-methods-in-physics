#ifndef NUMERICAL_METHODS_IN_PHYSICS_GAUSSSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_GAUSSSOLVER_H

#pragma once
#include "Base/IIntegralSolver.h"
#include <vector>
#include <string>

struct GaussTableEntry {
    std::vector<double> nodes;
    std::vector<double> weights;
};

static GaussTableEntry getGaussTable(size_t n) {
    switch (n) {
        case 2:
            return {{-0.5773502691896257, 0.5773502691896257},
                    { 1.0000000000000000, 1.0000000000000000}};
        case 3:
            return {{-0.7745966692414834, 0.0, 0.7745966692414834},
                    { 0.5555555555555556,  0.8888888888888888,
                        0.5555555555555556}};
        case 4:
            return {{-0.8611363115940526, -0.3399810435848563,
                      0.3399810435848563,  0.8611363115940526},
                    { 0.3478548451374539,  0.6521451548625461,
                      0.6521451548625461,  0.3478548451374539}};
        case 5:
            return {{-0.9061798459386640, -0.5384693101056831, 0.0,
                      0.5384693101056831,  0.9061798459386640},
                    { 0.2369268850561891, 0.4786286704993665,
                        0.5688888888888889,
                        0.4786286704993665,  0.2369268850561891}};
        case 6:
            return {{-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
                      0.2386191860831969,  0.6612093864662645,  0.9324695142031521},
                    { 0.1713244923791704,  0.3607615730481386,  0.4679139345726910,
                      0.4679139345726910,  0.3607615730481386,  0.1713244923791704}};
        case 7:
            return {{-0.9491079123427585, -0.7415311855993945, -0.4058451513773972,
                      0.0,
                      0.4058451513773972,  0.7415311855993945,  0.9491079123427585},
                    { 0.1294849661688697,  0.2797053914892766,  0.3818300505051189,
                      0.4179591836734694,
                      0.3818300505051189,  0.2797053914892766,  0.1294849661688697}};
        case 8:
            return {{-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                      -0.1834346424956498,
                       0.1834346424956498,  0.5255324099163290,  0.7966664774136267,
                       0.9602898564975363},
                    { 0.1012285362903763,  0.2223810344533745,  0.3137066458778873,
                      0.3626837833783620,
                      0.3626837833783620,  0.3137066458778873,  0.2223810344533745,
                      0.1012285362903763}};
        case 9:
            return {{-0.9681602395076261, -0.8360311073266358, -0.6133714327005904,
                      -0.3242534234038089,  0.0,
                       0.3242534234038089,  0.6133714327005904,  0.8360311073266358,
                       0.9681602395076261},
                    { 0.0812743883615744,  0.1806481606948574,  0.2606106964029354,
                      0.3123470770400029,  0.3302393550012598,
                      0.3123470770400029,  0.2606106964029354,  0.1806481606948574,
                      0.0812743883615744}};
        case 10:
            return {{-0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
                      -0.4333953941292472, -0.1488743389816312,
                       0.1488743389816312,  0.4333953941292472,  0.6794095682990244,
                       0.8650633666889845,  0.9739065285171717},
                    { 0.0666713443086881,  0.1494513491505806,  0.2190863625159820,
                      0.2692667193099963,  0.2955242247147529,
                      0.2955242247147529,  0.2692667193099963,  0.2190863625159820,
                      0.1494513491505806,  0.0666713443086881}};
        default:
            // fallback: 2-point rule
            return {{-0.5773502691896257,  0.5773502691896257},
                    { 1.0000000000000000,  1.0000000000000000}};
    }
}

class GaussSolver : public IIntegralSolver {
public:
    GaussSolver(size_t max_nodes = 10) : max_nodes_(std::max<size_t>(2, std::min<size_t>(10, max_nodes))) {}
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