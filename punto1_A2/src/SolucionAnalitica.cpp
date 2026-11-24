#include "SolucionAnalitica.h"

#include <cmath>

namespace {
    constexpr double kPi = 3.14159265358979323846;
}

SolucionAnalitica::SolucionAnalitica(double x0, double sigma, double k0, double hbar, double masa):
                    x0_(x0), sigma_(sigma), k0_(k0), hbar_(hbar), masa_(masa) {}

std::complex<double> SolucionAnalitica::Libre(double x, double t, double xc,
                                              double kc) const {
  const std::complex<double> unidad_imaginaria(0.0, 1.0);
  const double tau = hbar_ * t / (2.0 * masa_ * sigma_ * sigma_);
  const std::complex<double> factor = 1.0 + unidad_imaginaria * tau;
  const double desplazamiento = x - xc - hbar_ * kc * t / masa_;
  const double normalizacion = std::pow(2.0 * kPi * sigma_ * sigma_, -0.25);
  const std::complex<double> fase = unidad_imaginaria *
      (kc * x - hbar_ * kc * kc * t / (2.0 * masa_));
  return normalizacion / std::sqrt(factor) *
         std::exp(fase - desplazamiento * desplazamiento /
                             (4.0 * sigma_ * sigma_ * factor));
}

std::complex<double> SolucionAnalitica::Evaluar(double x, double t) const {

    return Libre(x, t, x0_, k0_) - Libre(x, t, -x0_, -k0_);
}