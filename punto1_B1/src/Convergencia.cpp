#include "Convergencia.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "CrankNicolson.h"
#include "EstadoCuantico.h"
#include "Hamiltoniano.h"
#include "Observables.h"
#include "PotencialArmonico.h"
#include "SolucionAnalitica.h"

namespace {

double ErrorL2(const Eigen::VectorXd& a, const Eigen::VectorXd& b, double dx) {
  return std::sqrt((a - b).squaredNorm() * dx);
}
}  // namespace

Convergencia::Convergencia(double x_min, double x_max, double masa,
                           double omega, double x0, double t_comparacion)
    : x_min_(x_min),
      x_max_(x_max),
      masa_(masa),
      omega_(omega),
      x0_(x0),
      t_comparacion_(t_comparacion) {
  std::filesystem::create_directories("results/convergencia");
  std::filesystem::create_directories("results/graficas");
  std::filesystem::create_directories("scripts");
}

Convergencia::Resultado Convergencia::Ejecutar(int n_interior,
                                               double dt) const {
  Malla malla(x_min_, x_max_, n_interior);

  // Ancho del estado fundamental (hbar = 1): estado coherente.
  const double sigma = std::sqrt(0.5 / (masa_ * omega_));
  EstadoCuantico estado(malla);
  estado.InitGaussiana(x0_, sigma, 0.0);

  PotencialArmonico potencial(masa_, omega_);
  Hamiltoniano hamiltoniano(malla, potencial.Vector(malla), 1.0, masa_);
  CrankNicolson propagador(hamiltoniano, dt);
  Observables observables(hamiltoniano);

  Resultado r;
  r.dx = malla.dx();
  r.x = malla.nodos();
  r.energia_inicial = observables.EnergiaMedia(estado);

  const int n_pasos = static_cast<int>(std::round(t_comparacion_ / dt));
  for (int paso = 0; paso < n_pasos; ++paso) propagador.Avanzar(estado);

  r.rho = observables.Densidad(estado);
  r.energia_final = observables.EnergiaMedia(estado);
  r.norma_final = estado.Norma();
  return r;
}

void Convergencia::EstudioEspacial(const std::vector<int>& n_interiores,
                                   double dt) {
  const SolucionAnalitica analitica(masa_, omega_, x0_);

  std::ofstream salida("results/convergencia/espacial.dat");
  salida << std::scientific << std::setprecision(8);
  salida << "# dx  error_L2_rho  |<H>-<H>_exacta|  |dH|_deriva  |N-1|\n";

  for (int n : n_interiores) {
    const Resultado r = Ejecutar(n, dt);

    Eigen::VectorXd rho_exacta(r.x.size());
    for (int j = 0; j < r.x.size(); ++j) {
      rho_exacta(j) = analitica.Densidad(r.x(j), t_comparacion_);
    }

    salida << r.dx << " " << ErrorL2(r.rho, rho_exacta, r.dx) << " "
           << std::abs(r.energia_inicial - analitica.EnergiaMedia()) << " "
           << std::abs(r.energia_final - r.energia_inicial) << " "
           << std::abs(r.norma_final - 1.0) << "\n";
    std::cout << "Espacial: n = " << n << ", dx = " << r.dx << " listo\n";
  }
}

void Convergencia::EstudioTemporal(int n_interior,
                                   const std::vector<double>& dts) {
    const Resultado referencia = Ejecutar(n_interior, dts.back() / 8.0);

  std::ofstream salida("results/convergencia/temporal.dat");
  salida << std::scientific << std::setprecision(8);
  salida << "# dt  error_L2_rho(vs dt_ref)  |dH|_deriva  |N-1|\n";

  for (double dt : dts) {
    const Resultado r = Ejecutar(n_interior, dt);
    salida << dt << " " << ErrorL2(r.rho, referencia.rho, r.dx) << " "
           << std::abs(r.energia_final - r.energia_inicial) << " "
           << std::abs(r.norma_final - 1.0) << "\n";
    std::cout << "Temporal: dt = " << dt << " listo\n";
  }
}

void Convergencia::Graficar() const {
  std::ofstream gp("scripts/graficar_convergencia.plt");
  gp << "# Generado por Convergencia::Graficar()\n";
  gp << "T = " << t_comparacion_ << "\n";
  gp << R"GP(
set terminal pngcairo size 900,600 enhanced
set grid
set logscale xy
set format xy '%g'
set fit quiet
set fit logfile 'results/convergencia/fit.log'

# Ajuste lineal en log-log: log10(error) = p * log10(h) + b
f_esp(x) = p_esp*x + b_esp
f_ene(x) = p_ene*x + b_ene
f_tem(x) = p_tem*x + b_tem
fit f_esp(x) 'results/convergencia/espacial.dat' using (log10($1)):(log10($2)) via p_esp,b_esp
fit f_ene(x) 'results/convergencia/espacial.dat' using (log10($1)):(log10($3)) via p_ene,b_ene
fit f_tem(x) 'results/convergencia/temporal.dat' using (log10($1)):(log10($2)) via p_tem,b_tem

set output 'results/graficas/convergencia_espacial.png'
set xlabel '{/Symbol D}x'
set ylabel 'error L2 de {/Symbol r}(x,T) frente a la solucion analitica'
set title sprintf('Refinamiento espacial en t = %.1f', T)
plot 'results/convergencia/espacial.dat' using 1:2 with linespoints pointtype 7 linewidth 2 \
     title sprintf('numerico (pendiente = %.2f)', p_esp)

set output 'results/graficas/convergencia_energia.png'
set ylabel '|<H> - <H>_{exacta}|'
set title 'Error de la energia media frente a {/Symbol D}x'
plot 'results/convergencia/espacial.dat' using 1:3 with linespoints pointtype 7 linewidth 2 \
     title sprintf('numerico (pendiente = %.2f)', p_ene)

set output 'results/graficas/convergencia_temporal.png'
set xlabel '{/Symbol D}t'
set ylabel 'error L2 de {/Symbol r}(x,T) frente a la corrida de referencia'
set title sprintf('Refinamiento temporal en t = %.1f', T)
plot 'results/convergencia/temporal.dat' using 1:2 with linespoints pointtype 7 linewidth 2 \
     title sprintf('numerico (pendiente = %.2f)', p_tem)
)GP";
  gp.close();

  const int codigo = std::system("gnuplot scripts/graficar_convergencia.plt");
  if (codigo != 0) {
    std::cerr << "Aviso: no se pudo ejecutar gnuplot (codigo " << codigo
              << "). Script en scripts/graficar_convergencia.plt\n";
  } else {
    std::cout << "Graficas de convergencia en results/graficas/\n";
  }
}