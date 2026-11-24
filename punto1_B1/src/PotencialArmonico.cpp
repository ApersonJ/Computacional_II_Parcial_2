#include "PotencialArmonico.h"

PotencialArmonico::PotencialArmonico(double masa, double omega)
    : masa_(masa), omega_(omega) {}

double PotencialArmonico::Evaluar(double x) const {
  return 0.5 * masa_ * omega_ * omega_ * x * x;
}

Eigen::VectorXd PotencialArmonico::Vector(const Malla& malla) const {
  Eigen::VectorXd v(malla.n());
  for (int j = 0; j < malla.n(); ++j) v(j) = Evaluar(malla.x(j));
  return v;
}