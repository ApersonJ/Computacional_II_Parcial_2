#ifndef CRANKNICOLSON_H_
#define CRANKNICOLSON_H_

#include <complex>

#include "eigen.h"
#include "Hamiltoniano.h"
#include "EstadoCuantico.h"

class CrankNicolson {

    public:
        CrankNicolson(const Hamiltoniano& hamiltoniano, double dt);

        void Preparar();
        void Avanzar(EstadoCuantico& estado);

        double dt() const{ return dt_;}

    private:
        const Hamiltoniano& hamiltoniano_;
        double dt_;
        Eigen::SparseMatrix<std::complex<double>> A_;
        Eigen::SparseMatrix<std::complex<double>> B_;
        Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solucionador_;

};


#endif  // CRANKNICOLSON_H_
