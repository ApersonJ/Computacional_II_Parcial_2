#ifndef CONVERGENCIA_H_
#define CONVERGENCIA_H_

#include <vector>

#include "malla.h"

class Convergencia {
 public:
  Convergencia(double x_min, double x_max, double masa, double omega,
               double x0, double t_comparacion);

  void EstudioEspacial(const std::vector<int>& n_interiores, double dt);
  void EstudioTemporal(int n_interior, const std::vector<double>& dts);
  void Graficar() const;

 private:
  struct Resultado {
    double dx;
    Eigen::VectorXd x;
    Eigen::VectorXd rho;
    double energia_inicial;
    double energia_final;
    double norma_final;
  };

  Resultado Ejecutar(int n_interior, double dt) const;

  double x_min_;
  double x_max_;
  double masa_;
  double omega_;
  double x0_;
  double t_comparacion_;
};

#endif  // CONVERGENCIA_H_

