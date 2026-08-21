#ifndef NUMERICAL_METHODS_IN_PHYSICS_SCALARFIELD_H
#define NUMERICAL_METHODS_IN_PHYSICS_SCALARFIELD_H
#pragma once
#include <vector>
#include <cassert>
#include <algorithm>

class ScalarField {
public:
    ScalarField() = default;
    ScalarField(int nx, int ny, double value = 0.0) { resize(nx, ny, value); }

    void resize(int nx, int ny, double value = 0.0) {
        nx_ = nx;
        ny_ = ny;
        data_.assign(static_cast<size_t>(nx_) * static_cast<size_t>(ny_), value);
    }

    int nx() const { return nx_; }
    int ny() const { return ny_; }

    void fill(double value) {
        std::fill(data_.begin(), data_.end(), value);
    }

    double& operator()(int i, int j) {
        assert(i >= 0 && i < nx_);
        assert(j >= 0 && j < ny_);
        return data_[static_cast<size_t>(j) * nx_ + i];
    }

    double operator()(int i, int j) const {
        assert(i >= 0 && i < nx_);
        assert(j >= 0 && j < ny_);
        return data_[static_cast<size_t>(j) * nx_ + i];
    }

private:
    int nx_ = 0;
    int ny_ = 0;
    std::vector<double> data_;
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_SCALARFIELD_H
