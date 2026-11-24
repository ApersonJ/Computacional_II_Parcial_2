#include "Simulacion.h"

#include <iomanip>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

Simulacion::Simulacion(const Malla& malla, EstadoCuantico& estado, CrankNicolson& propagador,
                const Observables& observables, double t_final, int pasos_serie, int pasos_foto) :

            malla_(malla), estado_(estado),
            propagador_(propagador), observables_(observables), t_final_(t_final),
            pasos_serie_(pasos_serie), pasos_foto_(pasos_foto), n_fotos_(0) {}

void Simulacion::Ejecutar() {
  std::filesystem::create_directories("results/graficas");
  std::filesystem::create_directories("scripts");

     std::ofstream serie("results/series.dat");
     std::ofstream densidad("results/densidad.dat");
     serie << std::setprecision(12);
     densidad << std::setprecision(8);

     serie << "# t  norma <x> <p> <H> prob_cerca_pared\n";
     densidad << "# x t rho\n";

     const double dt = propagador_.dt();
     const int n_pasos = static_cast<int>(std::round(t_final_ / dt));
     
     std::vector<double> t_serie, x_serie, px_serie;

    n_fotos_ = 0;
    for (int paso = 0; paso <= n_pasos; ++paso) {
        const double t = paso * dt;
        if (paso % pasos_serie_ == 0) {
            GuardarSerie(serie, t);
            t_serie.push_back(t);
            x_serie.push_back(observables_.PosicionMedia(estado_));
            px_serie.push_back(observables_.MomentumMedio(estado_));
        }
        if (paso % pasos_foto_ == 0) {
        GuardarFoto(densidad, t);
        ++n_fotos_;
    }
        if (paso < n_pasos) propagador_.Avanzar(estado_);
  }

  serie.close();     // Volcar a disco antes de que Gnuplot lea los datos.
  densidad.close();

  GuardarEhrenfest(t_serie, x_serie, px_serie);
  EscribirScriptGnuplot();
}

void Simulacion::GuardarSerie(std::ofstream& archivo, double t) const {

    archivo << t << " " << estado_.Norma() << " "
            << observables_.PosicionMedia(estado_) << " "
            << observables_.MomentumMedio(estado_) << " "
            << observables_.EnergiaMedia(estado_) << " "
            << observables_.ProbabilidadRegion(estado_, -5.0, 0.0) << "\n";
}

void Simulacion::GuardarFoto(std::ofstream& archivo, double t) const {
    const Eigen::VectorXd rho = observables_.Densidad(estado_);
    for (int j = 0; j < malla_.n(); ++j){
        archivo << malla_.x(j) << " " << t << " " << rho(j) << "\n";
    }
    archivo << "\n";
}

void Simulacion::GuardarEhrenfest(const std::vector<double>& t,
                              const std::vector<double>& x, const std::vector<double>& px) const {

 std::ofstream archivo("results/ehrenfest.dat");
 archivo << std::setprecision(12);
 archivo << "# t <x> <p> m(d<x>/dt) diferencia\n";

 const double masa = 1.0;
 double error_max = 0.0;
 
 for (std::size_t k = 1; k + 1 < t.size(); ++k){
  const double dxdt = masa*(x[k+1] - x[k - 1])/ (t[k + 1] - t[k - 1]);
  const double diferencia = px[k] - dxdt;
  archivo << t[k] << " " << px[k] << " " << dxdt << " " << diferencia << "\n";
  error_max = std::max(error_max, std::abs(diferencia));
 }

 std::cout << "Ehrenfest: max |<px> - m d<x>dt| = "<< error_max << "\n";
}

void Simulacion::EscribirScriptGnuplot() const {
  std::filesystem::create_directories("results/graficas");
  std::filesystem::create_directories("scripts");

  std::ofstream gp("scripts/graficar.plt");
  const double dt_foto = pasos_foto_ * propagador_.dt();

  gp << "# Generado por Simulacion::EscribirScriptGnuplot()\n";
  gp << "dt_foto = " << dt_foto << "\n";
  gp << "n_fotos = " << n_fotos_ << "\n";
  gp << R"GP(
set terminal pngcairo size 900,600 enhanced
set grid
set key top right

# --- Series temporales ---
set xlabel 't'

# Norma: casi constante en 1, rango muy estrecho y con decimales visibles.
set output 'results/graficas/norma.png'
set ylabel 'N(t)'
set yrange [1-1e-6:1+1e-6]
set format y '%.8f'
plot 'results/series.dat' using 1:2 with lines lw 2 title 'norma'
set format y '% h'
set autoscale y

set output 'results/graficas/x_medio.png'
set ylabel '<x>(t)'
plot 'results/series.dat' using 1:3 with lines lw 2 title '<x>'

set output 'results/graficas/px_medio.png'
set ylabel '<p_x>(t)'
plot 'results/series.dat' using 1:4 with lines lw 2 title '<p_x>'

# Energía: casi constante, se deja el rango automático.
set output 'results/graficas/energia.png'
set ylabel '<H>(t)'
set format y '%.6f'
plot 'results/series.dat' using 1:5 with lines lw 2 title '<H>'
set format y '% h'

set output 'results/graficas/prob_cerca_pared.png'
set ylabel 'P(-5 < x < 0)'
plot 'results/series.dat' using 1:6 with lines lw 2 title 'prob. cerca de la pared'

# --- Perfiles de densidad en varios instantes ---
set output 'results/graficas/densidad_perfiles.png'
set xlabel 'x'
set ylabel '{/Symbol r}(x,t) = |{/Symbol y}|^2'
plot for [i=0:n_fotos-1:4] 'results/densidad.dat' every :::i::i using 1:3 with lines lw 2 title sprintf('t = %.2f', i*dt_foto)

# --- Ehrenfest: <px> frente a m d<x>/dt ---
set output 'results/graficas/ehrenfest.png'
set xlabel 't'
set ylabel '<p_x>  y  m d<x>/dt'
plot 'results/ehrenfest.dat' using 1:2 with lines lw 3 title '<p_x> (derivada espacial)', \
     'results/ehrenfest.dat' using 1:3 with points pt 6 ps 0.8 title 'm d<x>/dt (derivada temporal)'

set output 'results/graficas/ehrenfest_diferencia.png'
set ylabel '<p_x> - m d<x>/dt'
plot 'results/ehrenfest.dat' using 1:4 with lines lw 2 title 'diferencia'

# --- Mapa de calor de rho(x,t) ---
set output 'results/graficas/mapa_densidad.png'
set xlabel 'x'
set ylabel 't'
set pm3d map
set palette defined (0 'white', 1 'blue', 2 'red')
splot 'results/densidad.dat' using 1:2:3 notitle
)GP";

  gp.close();  // El .plt debe estar completo en disco antes de ejecutarlo.

  const int codigo = std::system("gnuplot scripts/graficar.plt");
  if (codigo != 0) {
    std::cerr << "Aviso: no se pudo ejecutar gnuplot (codigo " << codigo
              << "). El script quedo en scripts/graficar.plt; "
              << "verifica que Gnuplot este instalado y en el PATH.\n";
  } else {
    std::cout << "Graficas generadas en results/graficas/\n";
  }
}


