#ifndef ESTADO_CUANTICO_H_
#define ESTADO_CUANTICO_H_

#include "eigen.h"
#include "malla.h"

class EstadoCuantico {

    public:
        EstadoCuantico(const Malla& malla);

        void InitGaussiana(double x0, double sigma, double k0);

        void Normalizar();

        double Norma() const;

        const Eigen::VectorXcd& amplitudes() const{
            return psi_;
        }
        Eigen::VectorXcd& amplitudes(){
            return psi_;
        }
        const Malla& malla() const{
            return malla_;
        }

    private:
        const Malla& malla_;
        Eigen::VectorXcd psi_;
};

#endif  // ESTADO_CUANTICO_H_