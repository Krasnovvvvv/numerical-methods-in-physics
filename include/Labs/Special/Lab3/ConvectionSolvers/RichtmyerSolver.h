#ifndef NUMERICAL_METHODS_IN_PHYSICS_RICHTMYERSOLVER_H
#define NUMERICAL_METHODS_IN_PHYSICS_RICHTMYERSOLVER_H

#pragma once

#include "Labs/Special/Lab3/Base/IConvectionSolver.h"

#include <vector>
#include <string>
#include <cstddef>

class RichtmyerSolver : public IConvectionSolver {
public:
    [[nodiscard]] std::string name() const override {
        return "Richtmyer (LW two-step, lecture-form)";
    }

    // ConvectionTask:
    //  - задаёт граничные значения u[0], u[N-1] до вызова step,
    //  - после step перезапишет u_next[0], u_next[N-1] по граничным условиям.
    //
    // Внутри:
    //   - нечётные вызовы step: выполняем 1-й шаг (Lax) j -> j+1
    //   - чётные вызовы step:  выполняем 2-й шаг (leapfrog) j,j+1 -> j+2
    //
    // Формулы:
    //   1-й шаг:
    //      U_{i,j+1} = 0.5 (U_{i+1,j} + U_{i-1,j})
    //                 - (c/2) (U_{i+1,j} - U_{i-1,j})
    //
    //   2-й шаг:
    //      U_{i,j+2} = U_{i,j} - c ( U_{i+1,j+1} - U_{i-1,j+1} )
    //
    void step(std::vector<double>&       u_next,
              const std::vector<double>& u,
              double                     c) const override
    {
        const std::size_t N = u.size();
        if (u_next.size() != N)
            u_next.resize(N);

        if (N < 3) {
            u_next = u;
            return;
        }

        if (!has_prev_) {
            // --- 1-й шаг: j -> j+1 (Lax) ---
            // prev_ пока пустой; считаем по текущему u (слой j)
            for (std::size_t i = 1; i + 1 < N; ++i) {
                u_next[i] =
                    0.5 * (u[i + 1] + u[i - 1])
                    - 0.5 * c * (u[i + 1] - u[i - 1]);
            }

            // границы не трогаем: ConvectionTask выставит их после вызова
            u_next[0]     = u[0];
            u_next[N - 1] = u[N - 1];

            // запоминаем слой j в prev_ и отмечаем, что следующий шаг будет leapfrog
            prev_      = u;
            has_prev_  = true;
            do_leap_   = true;
        } else if (do_leap_) {
            // --- 2-й шаг: j, j+1 -> j+2 (leapfrog) ---
            // prev_ = U_{i,j}, u = U_{i,j+1}
            const std::vector<double>& Uj   = prev_;
            const std::vector<double>& Ujp1 = u;

            for (std::size_t i = 1; i + 1 < N; ++i) {
                u_next[i] = Uj[i] - c * (Ujp1[i + 1] - Ujp1[i - 1]);
            }

            // границы снова оставляем на ConvectionTask
            u_next[0]     = u[0];
            u_next[N - 1] = u[N - 1];

            // теперь новый "j" для следующей пары шагов — это только что вычисленный слой j+2
            prev_     = u_next;
            do_leap_  = false;   // следующий вызов снова начнёт с Lax
        } else {
            // безопасный fallback: если по какой-то причине флаг сбился,
            // повторяем первый шаг (Lax) и восстанавливаем состояние.
            for (std::size_t i = 1; i + 1 < N; ++i) {
                u_next[i] =
                    0.5 * (u[i + 1] + u[i - 1])
                    - 0.5 * c * (u[i + 1] - u[i - 1]);
            }
            u_next[0]     = u[0];
            u_next[N - 1] = u[N - 1];

            prev_      = u;
            has_prev_  = true;
            do_leap_   = true;
        }
    }

private:
    // mutable, потому что step по интерфейсу логически const,
    // но нам нужно хранить внутреннее состояние метода
    mutable std::vector<double> prev_; // слой j
    mutable bool has_prev_ = false;    // есть ли сохранённый j
    mutable bool do_leap_  = false;    // следующий шаг — leapfrog (2-й шаг)
};

#endif // NUMERICAL_METHODS_IN_PHYSICS_RICHTMYERSOLVER_H
