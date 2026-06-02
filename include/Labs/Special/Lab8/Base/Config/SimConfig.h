#ifndef NUMERICAL_METHODS_IN_PHYSICS_SIMCONFIG_H
#define NUMERICAL_METHODS_IN_PHYSICS_SIMCONFIG_H
#pragma once

enum class InletType {
    Uniform,
    Parabolic
};

struct SimulationConfig {
    int Nx = 80;
    int Ny = 20;

    double Lx = 12.0;
    double Ly = 1.0;

    double rho = 1.0;
    double reynolds = 20.0;
    double uMean = 1.0;
    double pOutlet = 0.0;
    double vInlet = 0.0;
    InletType inletType = InletType::Uniform;

    int maxIterations = 1000;
    double alphaU = 0.7;
    double alphaV = 0.7;
    double alphaP = 1.0;

    double tolMass = 1e-4;
    double tolU = 1e-4;
    double tolV = 1e-4;
    double tolP = 1e-4;

    int printEvery = 25;

    double mu() const {
        return rho * uMean * Ly / reynolds;
    }
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_SIMCONFIG_H
