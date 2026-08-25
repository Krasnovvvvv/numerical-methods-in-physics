#ifndef NUMERICAL_METHODS_IN_PHYSICS_CASESETUP_H
#define NUMERICAL_METHODS_IN_PHYSICS_CASESETUP_H
#pragma once
#include "Labs/Special/Lab8/Base/Config/SimConfig.h"

inline SimulationConfig make_base_case() {
    SimulationConfig cfg;
    cfg.Nx = 80;
    cfg.Ny = 24;
    cfg.Lx = 12.0;
    cfg.Ly = 1.0;
    cfg.rho = 1.0;
    cfg.reynolds = 20.0;
    cfg.uMean = 1.0;
    cfg.pOutlet = 0.0;
    cfg.inletType = InletType::Uniform;
    cfg.maxIterations = 600;
    cfg.alphaU = 0.7;
    cfg.alphaV = 0.7;
    cfg.alphaP = 1.0;
    cfg.tolMass = 1e-4;
    cfg.tolU = 1e-4;
    cfg.tolV = 1e-4;
    cfg.tolP = 1e-4;
    cfg.printEvery = 20;
    return cfg;
}
#endif //NUMERICAL_METHODS_IN_PHYSICS_CASESETUP_H
