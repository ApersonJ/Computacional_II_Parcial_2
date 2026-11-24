#ifndef SOLUCIONANALITICA_H_
#define SOLUCIONANALITICA_H_

class SolucionAnalitica {
 public:
  SolucionAnalitica(double masa, double omega, double x0, double hbar = 1.0);

  // rho(x,t) = sqrt(m w / (pi hbar)) exp(-m w (x - x0 cos(w t))^2 / hbar)
  double Densidad(double x, double t) const;

  // <x>(t) = x0 cos(w t)
  double PosicionMedia(double t) const;

  // <H> = hbar w / 2 + m w^2 x0^2 / 2
  double EnergiaMedia() const;

 private:
  double masa_;
  double omega_;
  double x0_;
  double hbar_;
};

#endif  // SOLUCIONANALITICA_H_
