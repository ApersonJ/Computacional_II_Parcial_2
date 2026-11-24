#include <iostream>

#include "Convergencia.h"
#include "CrankNicolson.h"
#include "EstadoCuantico.h"
#include "Hamiltoniano.h"
#include "malla.h"
#include "Observables.h"
#include "Simulacion.h"


const double x_min = -25.0;   // frontera artificial izquierda
const double x_max = 0.0;     // pared infinita
const double x0 = -10.0;      // centro inicial del paquete
const double sigma = 1.0;     // ancho del paquete
const double k0 = 5.0;        // numero de onda (hacia la pared)

int main() {
  int opcion;
  std::cout << "0 = simulacion de A2, 1 = estudio de convergencia: ";
  std::cin >> opcion;

  if (opcion == 0) {
    Malla malla(x_min, x_max, 499);

    EstadoCuantico estado(malla);
    estado.InitGaussiana(x0, sigma, k0);

    Eigen::VectorXd potencial = Eigen::VectorXd::Zero(malla.n());
    Hamiltoniano hamiltoniano(malla, potencial);

    CrankNicolson propagador(hamiltoniano, 0.002);
    Observables observables(hamiltoniano);

    Simulacion simulacion(malla, estado, propagador, observables, 4.0, 25, 125);
    simulacion.Ejecutar();
  } else {
    Convergencia convergencia(x_min, x_max, x0, sigma, k0, 2.0);
    convergencia.EstudioEspacial({249, 499, 999, 1999}, 0.002);
    convergencia.EstudioTemporal(999, {0.08, 0.04, 0.02, 0.01});
    convergencia.Graficar();
    convergencia.CompararDensidad(499, 0.002, {1.0, 2.0, 3.0});
  }
  return 0;
}