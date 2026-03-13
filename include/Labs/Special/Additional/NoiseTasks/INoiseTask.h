#ifndef NUMERICAL_METHODS_IN_PHYSICS_INOISETASK_H
#define NUMERICAL_METHODS_IN_PHYSICS_INOISETASK_H

#pragma once

#include <string>

class INoiseTask {
public:
    virtual ~INoiseTask() = default;
    virtual void run() = 0;
    virtual std::string name() const = 0;
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_INOISETASK_H