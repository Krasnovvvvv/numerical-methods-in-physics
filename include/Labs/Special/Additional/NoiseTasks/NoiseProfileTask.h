#ifndef NUMERICAL_METHODS_IN_PHYSICS_NOISEPROFILETASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_NOISEPROFILETASK_H

#pragma once

#include <vector>
#include "Labs/Special/Additional/NoiseTasks/INoiseTask.h"
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"
#include "Helpers/Plotter.h"

class NoiseProfileTask : public INoiseTask {
    double a_;
    double b_;
    double dt_;
    std::size_t N_;
    std::vector<double> tau_list_;
    Plotter& plotter_;
    std::vector<DichotomicProfile> profiles_;

public:
    NoiseProfileTask(
        double a,
        double b,
        double dt,
        std::size_t N,
        std::vector<double> tau_list,
        Plotter& plotter)
        : a_(a), b_(b), dt_(dt), N_(N), tau_list_(tau_list), plotter_(plotter)
    {}

    std::string name() const override {
        return "Noise profile for different correlation times";
    }

    void run() override {
        profiles_.clear();
        profiles_.reserve(tau_list_.size());

        for (double tau_c : tau_list_) {
            DichotomicNoise noise(a_, tau_c, dt_);
            profiles_.emplace_back(noise.generate(N_));
        }

        for (const auto& p : profiles_) {
            std::string label = "tau_c = " + std::to_string(p.tau_c);
            plotter_.plot(
                p.t,
                p.s,
                label,
                "t",
                "{/Symbol s}(t)",
                false
            );
        }
    }

    const std::vector<DichotomicProfile>& profiles() const {
        return profiles_;
    }

};

#endif //NUMERICAL_METHODS_IN_PHYSICS_NOISEPROFILETASK_H