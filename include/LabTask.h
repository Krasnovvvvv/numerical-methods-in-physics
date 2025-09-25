#ifndef NUMERICAL_METHODS_IN_PHYSICS_LABTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_LABTASK_H

#pragma once
#include <utility>

#include "IDataGenerator.h"
#include "ISolver.h"
#include "Timer.h"
#include "Plotter.h"

class LabTask {
public:
    using ExpectedCurveFunc = std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&)>;

    LabTask(IDataGenerator& g, ISolver& s, Plotter& p,
            std::optional<ExpectedCurveFunc> expected = std::nullopt)
            : generator(g), solver(s), plotter(p), expected(std::move(expected)) {}

    virtual void run(const std::vector<size_t>& sizes) = 0;
    virtual ~LabTask() = default;
protected:
    IDataGenerator& generator;
    ISolver& solver;
    Plotter& plotter;
    std::optional<ExpectedCurveFunc> expected;
};


#endif //NUMERICAL_METHODS_IN_PHYSICS_LABTASK_H
