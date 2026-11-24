#ifndef HAMILTONIANO_H_
#define HAMILTONIANO_H_

#include <complex>

#include "eigen.h"
#include "malla.h"

class Hamiltoniano {

    public:
        Hamiltoniano(const Malla& malla, const Eigen:: VectorXd& potencial,
                     double hbar = 1.0, double masa = 1.0);

        void Ensamblar();

        bool EsSimetrico(double tolerancia = 1e-12) const;

        const Eigen::SparseMatrix<double>& Matriz() const{
            return H_;
        }

        Eigen::SparseMatrix<std::complex<double>> MatrizCompleja() const;

        double hbar() const{ return hbar_; }
        double masa() const{ return masa_; }
        const Malla& malla() const{ return malla_; }

    private:
        const Malla& malla_;
        const Eigen::VectorXd potencial_;
        double hbar_;
        double masa_;
        Eigen::SparseMatrix<double> H_;
        };

#endif  // HAMILTONIANO_H_