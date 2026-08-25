#ifndef NUMERICAL_METHODS_IN_PHYSICS_LINEARSYSTEM_H
#define NUMERICAL_METHODS_IN_PHYSICS_LINEARSYSTEM_H
#pragma once
#include <vector>

struct LinearSystem {
    int n = 0;
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> vals;
    std::vector<double> rhs;

    void reset(int size) {
        n = size;
        rows.clear();
        cols.clear();
        vals.clear();
        rhs.assign(static_cast<size_t>(size), 0.0);
    }

    void add(int r, int c, double v) {
        rows.push_back(r);
        cols.push_back(c);
        vals.push_back(v);
    }

    void add_rhs(int r, double v) {
        rhs[static_cast<size_t>(r)] += v;
    }
};
#endif //NUMERICAL_METHODS_IN_PHYSICS_LINEARSYSTEM_H