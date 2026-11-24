#include "Observables.h"

#include <complex>

Observables::Observables(const Hamiltoniano& hamiltoniano) : 
    hamiltoniano_(hamiltoniano), malla_(hamiltoniano.malla()) {}

Eigen::VectorXd Observables::Densidad(const EstadoCuantico& estado) const{
    Eigen::VectorXd rho = estado.amplitudes().cwiseAbs2();
    return rho;
}

double Observables::PosicionMedia(const EstadoCuantico& estado) const {
    const Eigen::VectorXd rho = Densidad(estado);
    return malla_.nodos().dot(rho) * malla_.dx();
}

double Observables::MomentumMedio(const EstadoCuantico& estado) const {
    const Eigen::VectorXcd& psi = estado.amplitudes();
    const int n = malla_.n();
    const double dx = malla_.dx();
    const std::complex<double> menos_i_hbar(0.0, -hamiltoniano_.hbar());

    std::complex<double> suma(0.0, 0.0);
    for (int j = 1; j < n - 1; ++j){
        const std::complex<double> anterior = (j > 0) ? psi(j - 1) : 0.0;
        const std::complex<double> siguiente = (j < n - 1) ? psi(j + 1) : 0.0;
        const std::complex<double> derivada = (siguiente - anterior) / (2.0 * dx);
        suma += std::conj(psi(j)) * menos_i_hbar * derivada;
    }
    return (suma*dx).real();
}

double Observables::EnergiaMedia(const EstadoCuantico& estado) const {
    const Eigen::SparseMatrix<double>& H = hamiltoniano_.Matriz();
    const Eigen::VectorXd re = estado.amplitudes().real();
    const Eigen::VectorXd im = estado.amplitudes().imag();
    const Eigen::VectorXd H_re = H * re;
    const Eigen::VectorXd H_im = H * im;
    return (re.dot(H_re) + im.dot(H_im)) * malla_.dx();
}

double Observables::ProbabilidadRegion(const EstadoCuantico& estado, double x_a, double x_b) const {
    const Eigen::VectorXd rho = Densidad(estado);
    double suma = 0.0;
    for (int j = 0; j < malla_.n(); ++j){
        if (malla_.x(j) >= x_a && malla_.x(j) <= x_b){
            suma += rho(j);
        }
    }
    return suma*malla_.dx();

}
