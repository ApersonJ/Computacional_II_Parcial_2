#include "graphics.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

namespace {
	const string DAT_INDIV   = "../results/trayectoria_individual.dat";
	const string DAT_PROM    = "../results/promedios_ensamble.dat";
	const string DIR_ANIM    = "../results/animations/";
	const string DIR_SCRIPTS = "../scripts/";

	// ejecuta un script de gnuplot y avisa si algo falla
	bool RunGnuplot(const string& script) {
		string cmd = "gnuplot \"" + script + "\"";
		int ret = system(cmd.c_str());
		if (ret != 0) {
			cerr << "Error al ejecutar gnuplot (" << script
			     << ")" << endl;
			return false;
		}
		return true;
	}

	// Genera y ejecuta una animacion de la columna ycol contra la columna xcol.
	// La columna 1 de los .dat siempre es t (se usa para el titulo de cada cuadro).
	void AnimateXY(const string& datafile, const string& scriptName, const string& gifName,
	               int xcol, int ycol, const string& title,
	               const string& xlabel, const string& ylabel, const string& legend) {

		string scriptPath = DIR_SCRIPTS + scriptName;
		ofstream gp(scriptPath.c_str());
		if (!gp.is_open()) {
			cerr << "No se pudo crear el script " << scriptPath << endl;
			return;
		}

		gp << "datafile = '" << datafile << "'\n";
		gp << "set terminal gif animate delay 4 loop 0 size 900,500\n";
		gp << "set output '" << DIR_ANIM << gifName << "'\n";

		gp << "# tiempo (columna 1) para el titulo de cada cuadro\n";
		gp << "stats datafile using 1 nooutput\n";
		gp << "N    = STATS_records\n";
		gp << "tmin = STATS_min\n";
		gp << "tmax = STATS_max\n";

		gp << "# rangos fijos para que los ejes no cambien entre cuadros\n";
		gp << "stats datafile using " << xcol << ":" << ycol << " nooutput\n";
		gp << "xmin = STATS_min_x\n";
		gp << "xmax = STATS_max_x\n";
		gp << "ymin = STATS_min_y\n";
		gp << "ymax = STATS_max_y\n";
		gp << "padx = 0.05*(xmax - xmin)\n";
		gp << "pady = 0.05*(ymax - ymin)\n";
		gp << "if (padx == 0) padx = 1\n";
		gp << "if (pady == 0) pady = 1\n";
		gp << "\n";
		gp << "set xrange [xmin-padx:xmax+padx]\n";
		gp << "set yrange [ymin-pady:ymax+pady]\n";
		gp << "set xlabel '" << xlabel << "'\n";
		gp << "set ylabel '" << ylabel << "'\n";
		gp << "set grid\n";
		gp << "set key top left\n";
		gp << "\n";
		gp << "# maximo 150 cuadros\n";
		gp << "nframes = (N-1 > 150) ? 150 : N-1\n";
		gp << "\n";
		gp << "do for [k=0:nframes] {\n";
		gp << "    idx = (k*(N-1))/nframes\n";
		gp << "    ti = tmin + (tmax-tmin)*idx/(N-1)\n";
		gp << "    set title sprintf('" << title << "  (t = %.2f)', ti)\n";
		gp << "    plot datafile every ::0::idx using " << xcol << ":" << ycol
		   << " with lines lw 2 lc rgb '#1f77b4' title '" << legend << "', "
		   << "datafile every ::idx::idx using " << xcol << ":" << ycol
		   << " with points pt 7 ps 1.8 lc rgb '#d62728' title 'actual'\n";
		gp << "}\n";
		gp << "\n";
		gp << "unset output\n";

		gp.close();

		if (RunGnuplot(scriptPath))
			cout << "Animacion guardada en " << DIR_ANIM << gifName << endl;
	}
}

//----------------------------------------------------------
//   graficas con trayectoria individual (sin promediar)
//----------------------------------------------------------
void Graphics::Position() {
	AnimateXY(DAT_INDIV, "position.gp", "position.gif", 1, 2,
	          "Posicion de la particula", "t", "x(t)", "x(t)");
}

void Graphics::PhaseSpace() {
	AnimateXY(DAT_INDIV, "phasespace.gp", "phasespace.gif", 2, 4,
	          "Espacio de fases", "x", "p", "(x, p)");
}

void Graphics::Velocity() {
	AnimateXY(DAT_INDIV, "velocity.gp", "velocity.gif", 1, 3,
	          "Velocidad de la particula", "t", "v(t)", "v(t)");
}

void Graphics::Momentum() {
	AnimateXY(DAT_INDIV, "momentum.gp", "momentum.gif", 1, 4,
	          "Momento de la particula", "t", "p(t)", "p(t)");
}

void Graphics::TotalEnergy() {
	AnimateXY(DAT_INDIV, "totalenergy.gp", "totalenergy.gif", 1, 7,
	          "Energia total", "t", "E(t)", "E(t)");
}

//----------------------------------------------------------
//    			Graficas con datos promediados
//----------------------------------------------------------
void Graphics::MeanPosition() {
	AnimateXY(DAT_PROM, "meanposition.gp", "meanposition.gif", 1, 2,
	          "Posicion media del ensamble", "t", "<x>(t)", "<x>(t)");
}

void Graphics::MeanKinetic() {
	AnimateXY(DAT_PROM, "meankinetic.gp", "meankinetic.gif", 1, 5,
	          "Energia cinetica media", "t", "<K>(t)", "<K>(t)");
}

void Graphics::MeanPotential() {
	AnimateXY(DAT_PROM, "meanpotential.gp", "meanpotential.gif", 1, 6,
	          "Energia potencial media", "t", "<U>(t)", "<U>(t)");
}
