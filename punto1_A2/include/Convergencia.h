#ifndef CONVERGENCIA_H_
#define CONVERGENCIA_H_

#include <vector>
#include "malla.h"

class Convergencia {

    public:
        Convergencia(double x_min, double x_max, double x0, double sigma, double k0, double t_comparacion);

        void EstudioEspacial(const std::vector<int>& n_interiores, double dt);
        void EstudioTemporal(int n_interior, const std::vector<double>& dts);
        void Graficar() const;
        void CompararDensidad(int n_interior, double dt,
                        const std::vector<double>& tiempos);

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
    double x0_;
    double sigma_;
    double k0_;
    double t_comparacion_;
};

#endif // CONVERGENCIA_H_

