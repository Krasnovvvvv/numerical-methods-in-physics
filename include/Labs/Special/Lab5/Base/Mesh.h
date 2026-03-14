#ifndef NUMERICAL_METHODS_IN_PHYSICS_MESH_H
#define NUMERICAL_METHODS_IN_PHYSICS_MESH_H
#pragma once

#include <vector>
#include <stdexcept>

struct Node {
    double r;
    double z;
    bool is_dirichlet = false;
    double dirichlet_value = 0.0;
};

struct Element {
    std::size_t n[3];
};

class Mesh {
public:
    std::vector<Node> nodes_;
    std::vector<Element> elements_;

    std::vector<std::size_t> bottom_edge_nodes_;

    std::size_t nr_ = 0;
    std::size_t nz_ = 0;

    void build(double H, std::size_t nr, std::size_t nz) {
        if (nr < 2 || nz < 2)
            throw std::invalid_argument("Nr or Nz must be greater than 2");

        nr_ = nr;
        nz_ = nz;

        nodes_.clear();
        elements_.clear();
        bottom_edge_nodes_.clear();

        double dr = 1.0 / static_cast<double>(nr - 1);
        double dz = H / static_cast<double>(nz - 1);

        nodes_.reserve(nr_ * nz_);
        for (std::size_t j = 0; j < nz_; ++j) {
            double z = dz * static_cast<double>(j);
            for (std::size_t i = 0; i < nr_; ++i) {
                double r = dr * static_cast<double>(i);

                Node node;
                node.r = r;
                node.z = z;
                node.is_dirichlet = false;
                node.dirichlet_value = 0.0;

                if (i == nr_ - 1) {
                    node.is_dirichlet = true;
                    node.dirichlet_value = 0.0;
                }

                if (j == nz_ - 1 ) {
                    node.is_dirichlet = true;
                }

                nodes_.push_back(node);

                if (j == 0)
                    bottom_edge_nodes_.emplace_back(index(i, j));

            }
        }

        for (std::size_t j = 0; j < nz_ - 1; ++j) {
            for (std::size_t i = 0; i < nr_ - 1; ++i) {
                std::size_t n00 = index(i, j);
                std::size_t n10 = index(i + 1, j);
                std::size_t n01 = index(i, j + 1);
                std::size_t n11 = index(i + 1, j + 1);

                elements_.emplace_back(Element{{n00, n10, n11}});
                elements_.emplace_back(Element{{n00, n11, n01}});
            }
        }
    }

    std::size_t index(std::size_t i, std::size_t j) const {
        return j * nr_ + i;
    }
};

#endif //NUMERICAL_METHODS_IN_PHYSICS_MESH_H