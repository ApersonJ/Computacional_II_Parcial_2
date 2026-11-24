#include "CrankNicolson.h"

#include <stdexcept>

CrankNicolson::CrankNicolson(const Hamiltoniano& hamiltoniano, double dt):
    hamiltoniano_(hamiltoniano), dt_(dt) {
        Preparar();
    }

void CrankNicolson::Preparar() {
    const int n = hamiltoniano_.malla().n();
    const std::complex<double> factor(0.0, 0.5*dt_/hamiltoniano_.hbar());

    Eigen::SparseMatrix<std::complex<double>> H = hamiltoniano_.MatrizCompleja();
    Eigen::SparseMatrix<std::complex<double>> I(n, n);
    I.setIdentity();

    A_ = I + factor * H;
    B_ = I - factor * H;

    A_.makeCompressed();
    solucionador_.analyzePattern(A_);
    solucionador_.factorize(A_);
    if (solucionador_.info() != Eigen::Success) {
        throw std::runtime_error("CrankNicolson: factorization falló");

    }
}

void CrankNicolson::Avanzar(EstadoCuantico& estado) {
    Eigen::VectorXcd lado_derecho = B_ * estado.amplitudes();
    estado.amplitudes() = solucionador_.solve(lado_derecho);
    if (solucionador_.info() != Eigen::Success) {
        throw std::runtime_error("CrankNicolson: solución falló");
    }
}