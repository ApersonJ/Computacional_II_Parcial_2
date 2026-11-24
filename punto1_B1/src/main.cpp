#include <cmath>
#include <iostream>

#include "Convergencia.h"
#include "CrankNicolson.h"
#include "EstadoCuantico.h"
#include "Hamiltoniano.h"
#include "malla.h"
#include "Observables.h"
#include "PotencialArmonico.h"
#include "Simulacion.h"

// Parametros de B1 (unidades hbar = 1).
const double x_min = -10.0;
const double x_max = 10.0;
const double masa = 1.0;
const double omega = 1.0;
const double x0 = 3.0;        // desplazamiento inicial respecto al equilibrio
const double pi = 3.14159265358979323846;

int main() {
  int opcion;
  std::cout << "0 = simulacion de B1, 1 = estudio de convergencia: ";
  std::cin >> opcion;

  if (opcion == 0) {
    Malla malla(x_min, x_max, 399); 

    
    const double sigma = std::sqrt(0.5 / (masa * omega));
    EstadoCuantico estado(malla);
    estado.InitGaussiana(x0, sigma, 0.0);

    PotencialArmonico potencial(masa, omega);
    Hamiltoniano hamiltoniano(malla, potencial.Vector(malla), 1.0, masa);

    CrankNicolson propagador(hamiltoniano, 0.002);
    Observables observables(hamiltoniano);

    Simulacion simulacion(malla, estado, propagador, observables,
                          4.0 * pi / omega, 25, 125, x0, omega);
    simulacion.Ejecutar();
  } else {
    
    Convergencia convergencia(x_min, x_max, masa, omega, x0, 6.0);
   
    convergencia.EstudioEspacial({199, 399, 799, 1599}, 0.002);
    convergencia.EstudioTemporal(399, {0.08, 0.04, 0.02, 0.01});
    convergencia.Graficar();
  }
  return 0;
}