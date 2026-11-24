#ifndef SIMULACION_H
#define SIMULACION_H

#include <fstream>
#include <string>
#include <vector>

#include "CrankNicolson.h"
#include "malla.h"
#include "EstadoCuantico.h"
#include "Observables.h"

class Simulacion {
    public:
        Simulacion (const Malla& malla, EstadoCuantico& estado, CrankNicolson& propagador, 
            const Observables& observables, double t_final, int pasos_serie, int pasos_foto);

        void Ejecutar();


    private:
        void GuardarSerie(std::ofstream& archivo, double t) const;
        void GuardarFoto(std::ofstream& archivo, double t) const;
        void EscribirScriptGnuplot() const;
        void GuardarEhrenfest(const std::vector<double>& t,
                              const std::vector<double>& x, const std::vector<double>& px) const;

        const Malla& malla_;
        EstadoCuantico& estado_;
        CrankNicolson& propagador_;
        const Observables& observables_;
        double t_final_;
        int pasos_serie_;
        int pasos_foto_;
        int n_fotos_;
};

#endif // SIMULACION_H