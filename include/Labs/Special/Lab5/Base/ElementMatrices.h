#ifndef NUMERICAL_METHODS_IN_PHYSICS_ELEMENTMATRICES_H
#define NUMERICAL_METHODS_IN_PHYSICS_ELEMENTMATRICES_H
#pragma once

#include "Mesh.h"
#include "CylinderHeatProblem.h"
#include <Eigen/Dense>

inline void compute_element_matrices(const Mesh& mesh,
                                     const CylinderHeatProblem& prob,
                                     const Element& elem,
                                     Eigen::Matrix3d& Ke,
                                     Eigen::Vector3d& Fe)
{
    // Координаты вершин треугольника
    double r1 = mesh.nodes_[elem.n[0]].r;
    double z1 = mesh.nodes_[elem.n[0]].z;

    double r2 = mesh.nodes_[elem.n[1]].r;
    double z2 = mesh.nodes_[elem.n[1]].z;

    double r3 = mesh.nodes_[elem.n[2]].r;
    double z3 = mesh.nodes_[elem.n[2]].z;

    double detJ = (r2 - r1)*(z3 - z1) - (r3 - r1)*(z2 - z1);
    double area = 0.5 * std::abs(detJ);

    // B_i = [dphi_i/dr, dphi_i/dz]
    Eigen::Vector2d B1, B2, B3;

    B1 <<  (z2 - z3), (r3 - r2);
    B2 <<  (z3 - z1), (r1 - r3);
    B3 <<  (z1 - z2), (r2 - r1);

    B1 /= (2.0 * area);
    B2 /= (2.0 * area);
    B3 /= (2.0 * area);

    double r_mid = (r1 + r2 + r3) / 3.0;
    double z_mid = (z1 + z2 + z3) / 3.0;

    double f_mid = prob.rhs(r_mid, z_mid);

    // Ke_ij = r_mid * (B_j ⋅ B_i) * area
    Ke.setZero();
    Eigen::Vector2d Bs[3] = {B1, B2, B3};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Ke(i,j) = r_mid * Bs[j].dot(Bs[i]) * area;
        }
    }

    // F_i ≈ f_mid * r_mid * area/3
    Fe.setZero();
    double coeff = f_mid * r_mid * area / 3.0;
    Fe(0) = coeff;
    Fe(1) = coeff;
    Fe(2) = coeff;
}

#endif //NUMERICAL_METHODS_IN_PHYSICS_ELEMENTMATRICES_H