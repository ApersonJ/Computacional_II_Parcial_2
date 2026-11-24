#ifndef OBSERVABLES_H_
#define OBSERVABLES_H_

#include "eigen.h"
#include "malla.h"
#include "EstadoCuantico.h"
#include "Hamiltoniano.h"

class Observables {
    public:
        Observables(const Hamiltoniano& hamiltoniano);

        Eigen::VectorXd Densidad(const EstadoCuantico& estado) const;
        double PosicionMedia(const EstadoCuantico& estado) const;
        double MomentumMedio(const EstadoCuantico& estado) const;
        double EnergiaMedia(const EstadoCuantico& estado) const;
        double ProbabilidadRegion( const EstadoCuantico& estado, double x_a, double x_b) const;


    private:
        const Hamiltoniano& hamiltoniano_;
        const Malla& malla_;
};

#endif  // OBSERVABLES_H_