#ifndef NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONDIFFUSIONPROBLEM_H
#define NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONDIFFUSIONPROBLEM_H
#pragma once

#include <stdexcept>
#include "Mesh.h"
#include "BoundaryCondition.h"
#include "Labs/Special/Lab7/TVDLimiter/TVDLimiter.h"

enum class ConvectionScheme {
    Upwind,
    TVD
};

struct ConvectionDiffusionProblem {
    Mesh mesh;

    double rho   = 1.0;
    double u     = 0.0;
    double v     = 0.0;
    double Gamma = 1.0;

    double Su = 0.0;
    double Sp = 0.0;

    BoundaryCondition left   = BoundaryCondition::dirichlet(0.0);
    BoundaryCondition right  = BoundaryCondition::outlet();
    BoundaryCondition bottom = BoundaryCondition::dirichlet(0.0);
    BoundaryCondition top    = BoundaryCondition::outlet();

    ConvectionScheme scheme = ConvectionScheme::TVD;
    TVDLimiterType limiter  = TVDLimiterType::VanLeer;
    double swebyBeta        = 1.5;

    explicit ConvectionDiffusionProblem(const Mesh& m) : mesh(m) {}

    void validate() const {
        if (rho <= 0.0) {
            throw std::runtime_error("rho must be positive.");
        }
        if (Gamma < 0.0) {
            throw std::runtime_error("Gamma must be non-negative.");
        }
        if (limiter == TVDLimiterType::Sweby && (swebyBeta < 1.0 || swebyBeta > 2.0)) {
            throw std::runtime_error("Sweby beta must be in [1, 2].");
        }
    }

    double cellVolume() const { return mesh.dx * mesh.dy; }
    double areaEW() const { return mesh.dy; }
    double areaSN() const { return mesh.dx; }

    double Fw() const { return rho * u * areaEW(); }
    double Fe() const { return rho * u * areaEW(); }
    double Fs() const { return rho * v * areaSN(); }
    double Fn() const { return rho * v * areaSN(); }

    double Dw() const { return Gamma * areaEW() / mesh.dx; }
    double De() const { return Gamma * areaEW() / mesh.dx; }
    double Ds() const { return Gamma * areaSN() / mesh.dy; }
    double Dn() const { return Gamma * areaSN() / mesh.dy; }

    double Dwb() const { return Gamma * areaEW() / (0.5 * mesh.dx); }
    double Deb() const { return Gamma * areaEW() / (0.5 * mesh.dx); }
    double Dsb() const { return Gamma * areaSN() / (0.5 * mesh.dy); }
    double Dnb() const { return Gamma * areaSN() / (0.5 * mesh.dy); }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONDIFFUSIONPROBLEM_H