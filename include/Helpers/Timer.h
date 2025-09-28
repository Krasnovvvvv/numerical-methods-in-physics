#ifndef NUMERICAL_METHODS_IN_PHYSICS_TIMER_H
#define NUMERICAL_METHODS_IN_PHYSICS_TIMER_H

#pragma once
#include <chrono>

template <typename Units = std::chrono::milliseconds>
class Timer {
public:
    using clock = std::chrono::steady_clock;
    Timer() : start(clock::now()) {}
    void reset() { start = clock::now(); }
    auto elapsed() const { return std::chrono::duration_cast<Units>(clock::now() - start).count(); }
private:
    std::chrono::time_point<clock> start;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_TIMER_H
