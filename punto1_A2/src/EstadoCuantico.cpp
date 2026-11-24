#include "../include/EstadoCuantico.h"

#include <cmath>
#include <complex>

EstadoCuantico::EstadoCuantico(const Malla& malla): 
    malla_(malla), psi_(Eigen::VectorXcd::Zero(malla.n())){}

void EstadoCuantico::InitGaussiana(double x0, double sigma, double k0) {

    const std::complex<double> img(0.0, 1.0);

    for (int j = 0; j < malla_.n(); ++j) {
        const double x = malla_.x(j);
        const double envolvente = std::exp(-(x-x0)*(x-x0)/(4.0*sigma*sigma));
        psi_(j) = envolvente * std::exp(img*k0*x);

    }
    Normalizar();
}

double EstadoCuantico::Norma() const {
    return psi_.squaredNorm() * malla_.dx();
}

void EstadoCuantico::Normalizar(){
    psi_ /= std::sqrt(Norma());
}