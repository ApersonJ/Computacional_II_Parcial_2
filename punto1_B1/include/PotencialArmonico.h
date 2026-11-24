#ifndef POTENCIALARMONICO_H_
#define POTENCIALARMONICO_H_

#include "malla.h"

class PotencialArmonico {
 public:
  PotencialArmonico(double masa, double omega);

  // V en una posición.
  double Evaluar(double x) const;

  // V evaluado en todos los nodos de la malla (un valor por nodo).
  Eigen::VectorXd Vector(const Malla& malla) const;

  double masa() const { return masa_; }
  double omega() const { return omega_; }

 private:
  double masa_;
  double omega_;
};

#endif  // POTENCIALARMONICO_H_