
#ifndef NUMERICAL_METHODS_IN_PHYSICS_LABTASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_LABTASK_H

#pragma once
#include "IDataGenerator.h"
#include "ISolver.h"
#include "Timer.h"
#include "Plotter.h"
class LabTask {
public:
    LabTask(IDataGenerator* g, ISolver* s, Plotter* p)
            : generator(g), solver(s), plotter(p) {}
    virtual void run(size_t n) {
        auto A = generator->generateMatrix(n);
        auto b = generator->generateVector(n);
        Timer<> t;
        auto result = solver->solve(A, b);
        auto elapsed = t.elapsed();
        plotter->plot(std::vector<double>{static_cast<double>(elapsed)}, "Time, ms");
    }
    virtual ~LabTask() = default;
protected:
    IDataGenerator* generator;
    ISolver* solver;
    Plotter* plotter;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_LABTASK_H
