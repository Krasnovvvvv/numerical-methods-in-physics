#ifndef NUMERICAL_METHODS_IN_PHYSICS_FLOWRESULT_H
#define NUMERICAL_METHODS_IN_PHYSICS_FLOWRESULT_H
#pragma once
#include "Labs/Special/Lab8/Base/StaggeredGrid.h"
#include "Labs/Special/Lab8/Base/Config/SimConfig.h"
#include "Labs/Special/Lab8/Base/Field/ScalarField.h"
#include "Labs/Special/Lab8/Base/Field/UField.h"
#include "Labs/Special/Lab8/Base/Field/VField.h"
#include "Residuals.h"

struct FlowResult {
    StaggeredGrid grid;
    SimulationConfig cfg;
    ScalarField p;
    UField u;
    VField v;
    Residuals residuals;
    int iterations = 0;
    double maxYPlus = 0.0;
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_FLOWRESULT_H