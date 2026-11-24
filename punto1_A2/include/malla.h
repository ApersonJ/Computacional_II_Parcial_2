#ifndef MALLA_H_
#define MALLA_H_

#include "eigen.h"

class Malla {

    public:
        Malla(double x_min, double x_max, double n_interior);

        int n() const{
            return n_;
        }

        double dx() const{
            return dx_;
        }

        double x(int j) const{
            return nodos_(j);
        }

        const Eigen::VectorXd& nodos() const{
            return nodos_;
        }

    private:
        double x_min_;
        double x_max_;
        int n_;
        double dx_;
        Eigen::VectorXd nodos_;
};

#endif  // MALLA_H_