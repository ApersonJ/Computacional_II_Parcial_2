#include "../include/malla.h"

Malla::Malla(double x_min, double x_max, double n_interior) : x_min_(x_min), x_max_(x_max), n_(n_interior) {
    dx_ = (x_max_ - x_min_) / (n_ + 1);
    nodos_.resize(n_);
    for (int j = 0; j < n_; ++j) {
        nodos_(j) = x_min_ + (j + 1) * dx_;
    }
}