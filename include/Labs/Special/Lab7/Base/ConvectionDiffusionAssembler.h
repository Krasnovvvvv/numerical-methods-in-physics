#ifndef NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONDIFFUSIONASSEMBLER_H
#define NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONDIFFUSIONASSEMBLER_H
#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ConvectionDiffusionProblem.h"
#include "LinearSystem.h"
#include "ScalarField.h"
#include "Labs/Special/Lab7/TVDLimiter/TVDLimiter.h"

class ConvectionDiffusionAssembler {
public:
    static LinearSystem build(const ConvectionDiffusionProblem& pr, const ScalarField& phiLag) {
        pr.validate();
        validateCompatible(pr.mesh, phiLag.mesh);

        const Mesh& m = pr.mesh;
        LinearSystem sys(m);

        const double Fw = pr.Fw();
        const double Fe = pr.Fe();
        const double Fs = pr.Fs();
        const double Fn = pr.Fn();

        const double Dw = pr.Dw();
        const double De = pr.De();
        const double Ds = pr.Ds();
        const double Dn = pr.Dn();

        const double V = pr.cellVolume();

        for (int j = 0; j < m.Ny; ++j) {
            for (int i = 0; i < m.Nx; ++i) {
                const size_t p = static_cast<size_t>(m.id(i, j));

                double aW = 0.0, aE = 0.0, aS = 0.0, aN = 0.0, aP = 0.0;
                double rhs = pr.Su * V;

                if (m.hasWest(i)) {
                    aW = Dw + std::max(Fw, 0.0);
                } else {
                    const double y = m.yc(j);
                    if (pr.left.type == BoundaryType::Dirichlet) {
                        const double phiB = pr.left.value(y);

                        aP  += pr.Dwb();
                        rhs += pr.Dwb() * phiB;

                        if (Fw > 0.0) {
                            aP  += Fw;
                            rhs += Fw * phiB;
                        }
                    } else if (Fw > 0.0) {
                        throw std::runtime_error("West outlet with inflow is invalid.");
                    }
                }

                if (m.hasEast(i)) {
                    aE = De + std::max(-Fe, 0.0);
                } else {
                    const double y = m.yc(j);
                    if (pr.right.type == BoundaryType::Dirichlet) {
                        const double phiB = pr.right.value(y);

                        aP  += pr.Deb();
                        rhs += pr.Deb() * phiB;

                        if (Fe < 0.0) {
                            aP  += (-Fe);
                            rhs += (-Fe) * phiB;
                        }
                    } else if (Fe < 0.0) {
                        throw std::runtime_error("East outlet with inflow is invalid.");
                    }
                }

                if (m.hasSouth(j)) {
                    aS = Ds + std::max(Fs, 0.0);
                } else {
                    const double x = m.xc(i);
                    if (pr.bottom.type == BoundaryType::Dirichlet) {
                        const double phiB = pr.bottom.value(x);

                        aP  += pr.Dsb();
                        rhs += pr.Dsb() * phiB;

                        if (Fs > 0.0) {
                            aP  += Fs;
                            rhs += Fs * phiB;
                        }
                    } else if (Fs > 0.0) {
                        throw std::runtime_error("South outlet with inflow is invalid.");
                    }
                }

                if (m.hasNorth(j)) {
                    aN = Dn + std::max(-Fn, 0.0);
                } else {
                    const double x = m.xc(i);
                    if (pr.top.type == BoundaryType::Dirichlet) {
                        const double phiB = pr.top.value(x);

                        aP  += pr.Dnb();
                        rhs += pr.Dnb() * phiB;

                        if (Fn < 0.0) {
                            aP  += (-Fn);
                            rhs += (-Fn) * phiB;
                        }
                    } else if (Fn < 0.0) {
                        throw std::runtime_error("North outlet with inflow is invalid.");
                    }
                }

                if (pr.scheme == ConvectionScheme::TVD) {
                    rhs += Fe * eastFaceCorrection(pr, phiLag, i, j);
                    rhs -= Fw * westFaceCorrection(pr, phiLag, i, j);
                    rhs += Fn * northFaceCorrection(pr, phiLag, i, j);
                    rhs -= Fs * southFaceCorrection(pr, phiLag, i, j);
                }

                aP += aW + aE + aS + aN + (Fe - Fw) + (Fn - Fs) - pr.Sp * V;

                sys.aW[p] = aW;
                sys.aE[p] = aE;
                sys.aS[p] = aS;
                sys.aN[p] = aN;
                sys.aP[p] = aP;
                sys.b[p]  = rhs;
            }
        }

        return sys;
    }

private:
    static constexpr double eps = 1e-14;

    static void validateCompatible(const Mesh& a, const Mesh& b) {
        if (a.Nx != b.Nx || a.Ny != b.Ny || a.Lx != b.Lx || a.Ly != b.Ly) {
            throw std::runtime_error("Incompatible meshes.");
        }
    }

    static double limiter(const ConvectionDiffusionProblem& pr, double r) {
        return TVDLimiter::eval(pr.limiter, r, pr.swebyBeta);
    }

    static double safeRatio(double num, double den) {
        if (std::abs(den) < eps) return 0.0;
        return num / den;
    }

    static double ghostWest(const ConvectionDiffusionProblem& pr, const ScalarField& phi, int j) {
        const double bc = pr.left.value(pr.mesh.yc(j));
        const double in = phi(0, j);
        return (pr.left.type == BoundaryType::Dirichlet) ? (2.0 * bc - in) : in;
    }

    static double ghostEast(const ConvectionDiffusionProblem& pr, const ScalarField& phi, int j) {
        const int i = pr.mesh.Nx - 1;
        const double bc = pr.right.value(pr.mesh.yc(j));
        const double in = phi(i, j);
        return (pr.right.type == BoundaryType::Dirichlet) ? (2.0 * bc - in) : in;
    }

    static double ghostSouth(const ConvectionDiffusionProblem& pr, const ScalarField& phi, int i) {
        const double bc = pr.bottom.value(pr.mesh.xc(i));
        const double in = phi(i, 0);
        return (pr.bottom.type == BoundaryType::Dirichlet) ? (2.0 * bc - in) : in;
    }

    static double ghostNorth(const ConvectionDiffusionProblem& pr, const ScalarField& phi, int i) {
        const int j = pr.mesh.Ny - 1;
        const double bc = pr.top.value(pr.mesh.xc(i));
        const double in = phi(i, j);
        return (pr.top.type == BoundaryType::Dirichlet) ? (2.0 * bc - in) : in;
    }

    static double sample(const ConvectionDiffusionProblem& pr, const ScalarField& phi, int i, int j) {
        const Mesh& m = pr.mesh;
        if (m.inside(i, j)) return phi(i, j);
        if (i == -1 && 0 <= j && j < m.Ny) return ghostWest(pr, phi, j);
        if (i == m.Nx && 0 <= j && j < m.Ny) return ghostEast(pr, phi, j);
        if (j == -1 && 0 <= i && i < m.Nx) return ghostSouth(pr, phi, i);
        if (j == m.Ny && 0 <= i && i < m.Nx) return ghostNorth(pr, phi, i);
        throw std::runtime_error("Only one ghost layer is supported.");
    }

    static double eastFaceCorrection(const ConvectionDiffusionProblem& pr, const ScalarField& phi, int i, int j) {
        if (!pr.mesh.hasEast(i) || std::abs(pr.Fe()) < eps) return 0.0;

        double phiC, phiD, phiU;
        if (pr.Fe() > 0.0) {
            phiC = phi(i, j);
            phiD = phi(i + 1, j);
            phiU = sample(pr, phi, i - 1, j);
        } else {
            phiC = phi(i + 1, j);
            phiD = phi(i, j);
            phiU = sample(pr, phi, i + 2, j);
        }

        const double d = phiD - phiC;
        const double r = safeRatio(phiC - phiU, d);
        return 0.5 * limiter(pr, r) * d;
    }

    static double westFaceCorrection(const ConvectionDiffusionProblem& pr, const ScalarField& phi, int i, int j) {
        if (!pr.mesh.hasWest(i) || std::abs(pr.Fw()) < eps) return 0.0;

        double phiC, phiD, phiU;
        if (pr.Fw() > 0.0) {
            phiC = phi(i - 1, j);
            phiD = phi(i, j);
            phiU = sample(pr, phi, i - 2, j);
        } else {
            phiC = phi(i, j);
            phiD = phi(i - 1, j);
            phiU = sample(pr, phi, i + 1, j);
        }

        const double d = phiD - phiC;
        const double r = safeRatio(phiC - phiU, d);
        return 0.5 * limiter(pr, r) * d;
    }

    static double northFaceCorrection(const ConvectionDiffusionProblem& pr, const ScalarField& phi, int i, int j) {
        if (!pr.mesh.hasNorth(j) || std::abs(pr.Fn()) < eps) return 0.0;

        double phiC, phiD, phiU;
        if (pr.Fn() > 0.0) {
            phiC = phi(i, j);
            phiD = phi(i, j + 1);
            phiU = sample(pr, phi, i, j - 1);
        } else {
            phiC = phi(i, j + 1);
            phiD = phi(i, j);
            phiU = sample(pr, phi, i, j + 2);
        }

        const double d = phiD - phiC;
        const double r = safeRatio(phiC - phiU, d);
        return 0.5 * limiter(pr, r) * d;
    }

    static double southFaceCorrection(const ConvectionDiffusionProblem& pr, const ScalarField& phi, int i, int j) {
        if (!pr.mesh.hasSouth(j) || std::abs(pr.Fs()) < eps) return 0.0;

        double phiC, phiD, phiU;
        if (pr.Fs() > 0.0) {
            phiC = phi(i, j - 1);
            phiD = phi(i, j);
            phiU = sample(pr, phi, i, j - 2);
        } else {
            phiC = phi(i, j);
            phiD = phi(i, j - 1);
            phiU = sample(pr, phi, i, j + 1);
        }

        const double d = phiD - phiC;
        const double r = safeRatio(phiC - phiU, d);
        return 0.5 * limiter(pr, r) * d;
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_CONVECTIONDIFFUSIONASSEMBLER_H