#ifndef SOLUCIONANALITICA_H_
#define SOLUCIONANALITICA_H_

#include <complex>

class SolucionAnalitica {

    public:
        SolucionAnalitica(double x0, double sigma, double k0, double hbar = 1.0, double masa = 1.0);

        std::complex<double> Evaluar(double x, double t) const;

    private:

        std::complex<double> Libre(double x, double t, double xc, double kc) const;

        double x0_;
        double sigma_;
        double k0_;
        double hbar_;
        double masa_;

    
};

#endif // SOLUCIONANALITICA_H_
