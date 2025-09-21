#ifndef NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H
#define NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H

#pragma once
#include <matplot/matplot.h>
#include <vector>
class Plotter {
public:
    virtual void plot(const std::vector<double>& data, const std::string& label) {
        auto p = matplot::plot(data);
        p->display_name(label);
        matplot::legend();
        matplot::show();
    }
};


#endif //NUMERICAL_METHODS_IN_PHYSICS_PLOTTER_H
