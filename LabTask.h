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
    virtual void run(size_t n) = 0;
    virtual ~LabTask() = default;
protected:
    IDataGenerator* generator;
    ISolver* solver;
    Plotter* plotter;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_LABTASK_H
