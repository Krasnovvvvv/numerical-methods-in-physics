#ifndef NUMERICAL_METHODS_IN_PHYSICS_STEPDEPENDENCESTRATEGY_H
#define NUMERICAL_METHODS_IN_PHYSICS_STEPDEPENDENCESTRATEGY_H

#pragma once
#include "Base/IODESolver.h"
#include "Labs/Lab4/Tasks/IODETaskStrategy.h"
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <optional>

class StepDependenceStrategy : public IODETaskStrategy {
public:
    StepDependenceStrategy(IODESolver* solver, std::vector<std::string> comp_names = {}, std::string time_name = "t")
        : solver_(solver), componentNames(std::move(comp_names)), timeName(std::move(time_name)), referenceSolutions(std::nullopt) {}

    void setReferenceSolutions(const std::vector<std::function<double(double)>>& refs) {
        referenceSolutions = refs;
    }

    void run(const std::function<std::vector<double>(double, const std::vector<double>&)>& rhs,
             const std::vector<double>& y0,
             double t0, double tn,
             double /*h*/, double tol,
             unsigned short graphNumber,
             Plotter* plotter,
             std::vector<double> steps,
             std::vector<double> /*tolerances*/,
             std::vector<IODESolver*> /*solvers*/) override
    {
        if (!solver_ || !plotter) throw std::runtime_error("Solver or Plotter not set!");
        if (steps.empty()) throw std::invalid_argument("Steps vector is empty!");

        std::vector<std::vector<double>> all_xs, all_ys;
        std::vector<std::string> all_labels;

        for (double step : steps) {
            auto res = solver_->solve(rhs, y0, t0, tn, step, tol);
            size_t dim = res.y.empty() ? 0 : res.y[0].size();

            std::vector<std::vector<double>> y_components(dim);
            for (const auto& state : res.y)
                for (size_t j = 0; j < dim; ++j)
                    y_components[j].push_back(state[j]);

            if (graphNumber == 1) {
                // Solutions
                for (size_t j = 0; j < dim; ++j) {
                    std::string cname = (j < componentNames.size())
                        ? componentNames[j] : ("U" + std::to_string(j+1));
                    std::ostringstream oss;
                    oss << cname << ", h=" << step;
                    all_xs.push_back(res.t);
                    all_ys.push_back(y_components[j]);
                    all_labels.push_back(oss.str());
                }
            }

            if (graphNumber == 2) {
                // Derivations
                for (size_t j = 0; j < dim; ++j) {
                    std::vector<double> deriv, t_shifted;
                    for (size_t i = 1; i < res.t.size(); ++i) {
                        deriv.push_back((y_components[j][i] - y_components[j][i-1]) / (res.t[i] - res.t[i-1]));
                        t_shifted.push_back(res.t[i]);
                    }
                    std::string cname = (j < componentNames.size())
                        ? componentNames[j] : ("U" + std::to_string(j+1));
                    std::ostringstream oss;
                    oss << "d/dt(" << cname << "), h=" << step;
                    all_xs.push_back(t_shifted);
                    all_ys.push_back(deriv);
                    all_labels.push_back(oss.str());
                }
            }

            if (graphNumber == 3) {
                // Euclidean error
                if (referenceSolutions.has_value() && !referenceSolutions->empty()) {
                    std::vector<double> euclid_errors(res.t.size(), NAN);
                    for (size_t i = 0; i < res.t.size(); ++i) {
                        double sum_sq = 0.0;
                        for (size_t comp = 0; comp < dim && comp < referenceSolutions->size(); ++comp) {
                            double val_num = res.y[i][comp];
                            double val_ex = (*referenceSolutions)[comp](res.t[i]);
                            double delta = val_num - val_ex;
                            sum_sq += delta * delta;
                        }
                        euclid_errors[i] = std::sqrt(sum_sq);
                    }
                    std::ostringstream oss;
                    oss << "Euclidean error, h=" << step;
                    all_xs.push_back(res.t);
                    all_ys.push_back(euclid_errors);
                    all_labels.push_back(oss.str());
                } else {
                    std::ostringstream oss;
                    oss << "Error estimate, h=" << step;
                    all_xs.push_back(res.t);
                    all_ys.push_back(res.errorEstimates);
                    all_labels.push_back(oss.str());
                }
            }

            if (graphNumber == 4) {
                // Phase trajectories
                for (size_t a = 0; a < dim; ++a) {
                    for (size_t b = a+1; b < dim; ++b) {
                        std::ostringstream oss;
                        oss << ((a < componentNames.size()) ? componentNames[a] : ("U" + std::to_string(a+1)))
                            << " vs "
                            << ((b < componentNames.size()) ? componentNames[b] : ("U" + std::to_string(b+1)))
                            << ", h=" << step;
                        all_xs.push_back(y_components[a]);
                        all_ys.push_back(y_components[b]);
                        all_labels.push_back(oss.str());
                    }
                }
            }
        }

        if (graphNumber == 1)
            plotter->plot(all_xs, all_ys, all_labels, timeName, "Ui");
        else if (graphNumber == 2)
            plotter->plot(all_xs, all_ys, all_labels, timeName, "dUi/dt");
        else if (graphNumber == 3)
            plotter->plot(all_xs, all_ys, all_labels, timeName, "Euclidean Error");
        else if (graphNumber == 4)
            plotter->plot(all_xs, all_ys, all_labels, "Ui", "Ui'");
    }

private:
    IODESolver* solver_;
    std::vector<std::string> componentNames;
    std::string timeName;
    std::optional<std::vector<std::function<double(double)>>> referenceSolutions;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_STEPDEPENDENCESTRATEGY_H