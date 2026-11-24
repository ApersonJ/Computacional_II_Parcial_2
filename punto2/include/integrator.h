#ifndef INTEGRATOR
#define INTEGRATOR

#include "particle.h"
#include "motormodel.h"
#include <random>

class Integrator {
private:
    std::mt19937 gen;
    std::normal_distribution<double> dist;
    
    MotorModel model;
    double tau_on, tau_off;

public:
	Integrator(unsigned int seed = 42) : gen(seed), dist(0.0, 1.0) {}
	
	void Inicie(MotorModel model_val, double tau_on_val, double tau_off_val){
		model = model_val;
		tau_on = tau_on_val;
		tau_off = tau_off_val;
	}
	
	//integrador estocástico, cuando gamma y kbT == 0 (caso determinista) pasa a ser Velocity-Verlet
	void GJFIntegrator(double t, double dt){
		double m     = model.getParticle().getM();
        double x_n   = model.getParticle().getX();
        double v_n   = model.getParticle().getV();
        double gamma = model.getGamma();
        double kBT   = model.getT(); // kb = 1 se asume en main
        
        //constantes del método GJF
        double b = 1/(1+gamma*dt/(2*m));
        double a = (1-gamma*dt/(2*m))/(1+gamma*dt/(2*m));
        
        //término de ruido gaussiano para la fuerza
        double xi = dist(gen); //numero aleatorio de una distribución gaussiana con media 0 y varianza 1
		double beta = sqrt(2*gamma*kBT*dt)*xi;
		
		//fuerza actual
		double F_n = model.force(t, tau_on, tau_off);
		
		//actualizar posición
		double x_next = x_n + b*dt*v_n + b*dt*dt/(2*m)*F_n + b*dt/(2*m)*beta;
		model.getParticle().setPosition(x_next);
		
		//actualizar fuerza
		double F_next = model.force(t+dt, tau_on, tau_off);
		
		//actualizar velocidad
		double v_next = a*v_n + dt/(2*m)*(a*F_n + F_next) + b/m*beta;
		model.getParticle().setVelocity(v_next);	
	}
	
	void setSeed(unsigned int newSeed) {
        gen.seed(newSeed);
        dist.reset();
    }
	
	MotorModel& getModel(){ return model;}
	
};

#endif
