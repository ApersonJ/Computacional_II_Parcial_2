#include "SolucionAnalitica.h"

#include <cmath>

namespace {
constexpr double kPi = 3.14159265358979323846;
}

SolucionAnalitica::SolucionAnalitica(double masa, double omega, double x0,
                                     double hbar)
    : masa_(masa), omega_(omega), x0_(x0), hbar_(hbar) {}

double SolucionAnalitica::Densidad(double x, double t) const {
  const double centro = PosicionMedia(t);
  const double a = masa_ * omega_ / hbar_;
  return std::sqrt(a / kPi) * std::exp(-a * (x - centro) * (x - centro));
}

double SolucionAnalitica::PosicionMedia(double t) const {
  return x0_ * std::cos(omega_ * t);
}

double SolucionAnalitica::EnergiaMedia() const {
  return 0.5 * hbar_ * omega_ + 0.5 * masa_ * omega_ * omega_ * x0_ * x0_;
}