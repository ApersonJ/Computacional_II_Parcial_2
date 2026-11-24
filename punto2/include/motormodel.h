#ifndef MOTORMODEL
#define MOTORMODEL

#include "particle.h"

class MotorModel {
private:
    Particle particle;
    double gamma;
    double T;
    
    double L, U0, alpha;

public:
	void Inicie(Particle particle_val, double gamma_val, double T_val, double L_val, double U0_val, double alpha_val){
		particle = particle_val;
		gamma = gamma_val;
		T = T_val;
		L = L_val;
		U0 = U0_val;
		alpha = alpha_val;
	}
	
	double force(double t, double tau_on, double tau_off);
	double potentialEnergy(double t, double tau_on, double tau_off);
	double kineticEnergy(){
		return particle.kineticEnergy();
	}

	double totalEnergy(double t, double tau_on, double tau_off) {
    	return kineticEnergy() + potentialEnergy(t, tau_on, tau_off);
	}

    double getGamma() const { return gamma; }
    double getT() const { return T; }
    Particle& getParticle() { return particle; }
};

#endif
