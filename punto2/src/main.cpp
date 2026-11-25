#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "particle.h"
#include "motormodel.h"
#include "integrator.h"
#include "graphics.h"

using namespace std;

//es para generar los bucles
void Simulation(int N_pasos, int N_ens, int saveStride, 
				double x0, double v0, double m, Integrator integrator, 
				double t0, double dt, double tau_on, double tau_off);

int main(){
	
	double x0, v0;									//condiciones iniciales
	double tau_on, tau_off; 								//tiempos de ensendido y apagado para el potencial diente de cierra
	double m, L, U0, alpha, gamma, T, kb=1;				//masa, periodo del potencial, profuncidad del poso, parï¿½metro de asimetrï¿½a, coeficiente de fricciï¿½n, temperatura
	double t0 = 0, tf, dt=0.001;							//tiempos
	int N_pasos, N_ens, saveStride = 10;	//discretización de la simulación, número de simulaciones para el ruido, cada cuanto guardar datos
	
	
	cout<<"Considerar friccion y ruido termico?"<<endl;
	cout<<"1. Si (modelo estocastico)"<<endl;
	cout<<"2. No (modelo determinista)"<<endl;
	int elec;
	cout<<"Eleccion: "; cin>>elec;
	
	while(elec != 1 && elec != 2){
		cout<<"Eleccion invalida, elija entre 1 o 2: "; cin>>elec;
	}
	
	cout<<"Duracion de la simulacion: "; cin>>tf;
	N_pasos = round((tf-t0)/dt);
	
	cout<<"Duracion de potencial encendido: "; cin>>tau_on;
	cout<<"Duracion de potencial apagado: "; cin>>tau_off;
	
	cout<<"Posicion inicial: "; cin>>x0;
	cout<<"Velocidad inicial: "; cin>>v0;
	cout<<"Masa: "; cin>>m;
	cout<<"Periodo del potencial: "; cin>>L;
	cout<<"Profundidad del poso: "; cin>>U0;
	cout<<"Parametro de asimetria (0, 0.5): "; cin>>alpha;
	
	if(elec == 1){
		cout<<"Coeficiente de friccion: "; cin>>gamma;
		cout<<"Energia termica: "; cin>>T;
		 N_ens = 500;
	}else{
		gamma = 0;
		T = 0;
		N_ens = 1;
	}
	
	//instanciar clases
	Particle particle;
	particle.Inicie(x0, v0, m);
	
	MotorModel model;
	model.Inicie(particle, gamma, T, L, U0, alpha);
	
	Integrator integrator;
	integrator.Inicie(model, tau_on, tau_off);
	
	Simulation(N_pasos, N_ens, saveStride, 
				x0, v0, m, integrator, 
				t0, dt, tau_on, tau_off);
	
	
	//realizar todas las animaciones			
	Graphics graphics;
	graphics.Position();
	graphics.PhaseSpace();
	graphics.Velocity();
	graphics.Momentum();
	graphics.TotalEnergy();
	graphics.MeanPosition();
	graphics.MeanKinetic();
	graphics.MeanPotential();
	
	return 0;
}



void Simulation(int N_pasos, int N_ens, int saveStride, 
				double x0, double v0, double m, Integrator integrator, 
				double t0, double dt, double tau_on, double tau_off){
	
    int numGuardados = (N_pasos / saveStride) + 1;
    vector<double> t_vec(numGuardados, 0.0);
    vector<double> x_prom(numGuardados, 0.0);
    vector<double> v_prom(numGuardados, 0.0);
    vector<double> K_prom(numGuardados, 0.0);
    vector<double> U_prom(numGuardados, 0.0);
    vector<double> E_prom(numGuardados, 0.0);

    double suma_v_media = 0.0;

    ofstream fileIndiv("../results/trayectoria_individual.dat");
    fileIndiv << "# t\tx\tv\tp\tK\tU\tE\n";
    fileIndiv << scientific << setprecision(8);

	// repetir la simulación con valores de ruido distinto (si se eligió método determinista esto solo se hace una vez)
    for (int k = 0; k < N_ens; ++k) {
        integrator.getModel().getParticle().setPosition(x0);
        integrator.getModel().getParticle().setVelocity(v0);
        integrator.setSeed(1000 + k);

        // Estado inicial en t0 (índice 0)
        double K0 = integrator.getModel().kineticEnergy();
        double U0_val = integrator.getModel().potentialEnergy(t0, tau_on, tau_off);
        double E0 = K0 + U0_val;

        if (k == 0) {
            t_vec[0] = t0;
            fileIndiv << t0 << "\t" << x0 << "\t" << v0 << "\t"
                      << integrator.getModel().getParticle().momentum() << "\t"
                      << K0 << "\t" << U0_val << "\t" << E0 << "\n";
        }
        x_prom[0] += x0;
        v_prom[0] += v0;
        K_prom[0] += K0;
        U_prom[0] += U0_val;
        E_prom[0] += E0;

		//bucle temporal, guarda datos cada saveStride bucles
        for (int n = 0; n < N_pasos; ++n) {
            double t_actual = t0 + n * dt;
            double t_sig    = t0 + (n + 1) * dt;

            integrator.GJFIntegrator(t_actual, dt);

            if ((n + 1) % saveStride == 0) {
                int idx = (n + 1) / saveStride;

                double x = integrator.getModel().getParticle().getX();
                double v = integrator.getModel().getParticle().getV();
                double p = integrator.getModel().getParticle().momentum();
                double K = integrator.getModel().kineticEnergy();
                double U = integrator.getModel().potentialEnergy(t_sig, tau_on, tau_off);
                double E = K + U;

                if (k == 0) {
                    t_vec[idx] = t_sig;
                    fileIndiv << t_sig << "\t" << x << "\t" << v << "\t"
                              << p << "\t" << K << "\t" << U << "\t" << E << "\n";
                }

                x_prom[idx] += x;
                v_prom[idx] += v;
                K_prom[idx] += K;
                U_prom[idx] += U;
                E_prom[idx] += E;
            }
        }

        double x_final = integrator.getModel().getParticle().getX();
        suma_v_media += (x_final - x0) / (N_pasos * dt);
    }
    fileIndiv.close();

	//almacenar todo en el .dat
    ofstream fileProm("../results/promedios_ensamble.dat");
    fileProm << "# t \t x_mean \t v_mean \t p_mean \t K_mean \t U_mean \t E_mean \n";
    fileProm << scientific << setprecision(8);

    for (int idx = 0; idx < numGuardados; ++idx) {
        double xm = x_prom[idx] / N_ens;
        double vm = v_prom[idx] / N_ens;
        fileProm << t_vec[idx] << "\t"
                 << xm << "\t"
                 << vm << "\t"
                 << m * vm << "\t"
                 << K_prom[idx] / N_ens << "\t"
                 << U_prom[idx] / N_ens << "\t"
                 << E_prom[idx] / N_ens << "\n";
    }
    fileProm.close();

    cout << "\nSimulacion finalizada." << endl;
    cout << "Velocidad media estacionaria <v>: " << (suma_v_media / N_ens) << endl;
	
}
