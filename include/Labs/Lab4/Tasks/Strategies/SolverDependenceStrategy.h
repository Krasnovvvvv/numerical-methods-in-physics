#ifndef NUMERICAL_METHODS_IN_PHYSICS_SOLVERDEPENDENCESTRATEGY_H
#define NUMERICAL_METHODS_IN_PHYSICS_SOLVERDEPENDENCESTRATEGY_H

#pragma once
#include "Base/IODESolver.h"
#include "Labs/Lab4/Tasks/IODETaskStrategy.h"
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>

class SolverDependenceStrategy : public IODETaskStrategy {
public:
    SolverDependenceStrategy(std::vector<std::string> comp_names = {}, std::string time_name = "t")
        : componentNames(std::move(comp_names)), timeName(std::move(time_name)) {}

    void run(const std::function<std::vector<double>(double, const std::vector<double>&)>& rhs,
             const std::vector<double>& y0,
             double t0, double tn,
             double h, double tol,
             unsigned short graphNumber,
             Plotter* plotter,
             std::vector<double> /*steps*/,
             std::vector<double> /*tolerances*/,
             std::vector<IODESolver*> solvers) override
    {
        if (solvers.empty()) throw std::invalid_argument("Solvers vector is empty!");
        if (!plotter) throw std::runtime_error("Plotter not set!");

        std::vector<std::vector<double>> all_xs, all_ys;
        std::vector<std::string> all_labels;

        for (auto* solver : solvers) {
            auto res = solver->solve(rhs, y0, t0, tn, h, tol);
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
                    oss << cname << ", " << solver->name();
                    all_xs.push_back(res.t);
                    all_ys.push_back(y_components[j]);
                    all_labels.push_back(oss.str());
                }
            }

            if (graphNumber == 2) {
                // Derivatives
                for (size_t j = 0; j < dim; ++j) {
                    std::vector<double> deriv, t_shifted;
                    for (size_t i = 1; i < res.t.size(); ++i) {
                        deriv.push_back((y_components[j][i] - y_components[j][i-1]) / (res.t[i] - res.t[i-1]));
                        t_shifted.push_back(res.t[i]);
                    }
                    std::string cname = (j < componentNames.size())
                        ? componentNames[j] : ("U" + std::to_string(j+1));
                    std::ostringstream oss;
                    oss << "d/dt(" << cname << "), " << solver->name();
                    all_xs.push_back(t_shifted);
                    all_ys.push_back(deriv);
                    all_labels.push_back(oss.str());
                }
            }

            if (graphNumber == 3) {
                // Error
                std::ostringstream oss;
                oss << "Error, " << solver->name();
                all_xs.push_back(res.t);
                all_ys.push_back(res.errorEstimates);
                all_labels.push_back(oss.str());
            }

            if (graphNumber == 4) {
                // Phase trajectories
                for (size_t a = 0; a < dim; ++a) {
                    for (size_t b = a + 1; b < dim; ++b) {
                        std::ostringstream oss;
                        oss << ((a < componentNames.size()) ? componentNames[a] : ("U" + std::to_string(a+1)))
                            << " vs "
                            << ((b < componentNames.size()) ? componentNames[b] : ("U" + std::to_string(b+1)))
                            << ", " << solver->name();
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
            plotter->plot(all_xs, all_ys, all_labels, timeName, "Error");
        else if (graphNumber == 4)
            plotter->plot(all_xs, all_ys, all_labels, "Ui", "Ui'");
    }

private:
    std::vector<std::string> componentNames;
    std::string timeName;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_SOLVERDEPENDENCESTRATEGY_H