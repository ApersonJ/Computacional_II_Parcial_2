#include "Hamiltoniano.h"

#include <stdexcept>
#include <vector>


Hamiltoniano::Hamiltoniano(const Malla& malla, const Eigen::VectorXd& potencial,
                           double hbar, double masa):

                malla_(malla), potencial_(potencial), hbar_(hbar), masa_(masa){

                    if (potencial_.size() != malla_.n()){
                        throw std::invalid_argument("Hamiltoniano: el potencial debe tener un valor por nodo");

                    }
                    Ensamblar();
                }


void Hamiltoniano::Ensamblar(){

    const int n = malla_.n();
    const double dx = malla_.dx();
    const double alfa = hbar_*hbar_/(2.0*masa_*dx*dx);

    std::vector<Eigen::Triplet<double>> ternas;
    ternas.reserve(3*n);

    for (int j = 0; j < n; ++j){
        ternas.emplace_back(j, j, 2.0*alfa + potencial_(j));
        if (j > 0) ternas.emplace_back(j, j-1, -alfa);
        if (j < n-1) ternas.emplace_back(j, j+1, -alfa);
    }
    H_.resize(n, n);
    H_.setFromTriplets(ternas.begin(), ternas.end());
}

bool Hamiltoniano::EsSimetrico(double tolerancia) const{
    Eigen::SparseMatrix<double> traspuesta = H_.transpose();
    Eigen::SparseMatrix<double> diferencia = H_ - traspuesta;
    return diferencia.norm() < tolerancia;
}

Eigen::SparseMatrix<std::complex<double>> Hamiltoniano::MatrizCompleja() const{
    return H_.cast<std::complex<double>>();
}