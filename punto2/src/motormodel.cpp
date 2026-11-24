#include "motormodel.h"
#include "particle.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>

using namespace std;

//----------------------------------------------------------
//					Definir el Potencial
//----------------------------------------------------------
double MotorModel::potentialEnergy(double t, double tau_on, double tau_off){
	//Coordenadas reducidas
	double tau = tau_on + tau_off;
	double x_red = particle.getX() - L*floor(particle.getX()/L);
	double t_red = t - tau*floor(t/tau);
	//modulador
	double mod;
	//modulador
	if(0<=t_red && t_red<tau_on){
		mod = 1;
		
	}else{
		mod = 0;
	}
	
	double U;
	if(x_red<alpha*L){
		U = U0/(alpha*L) * x_red;
	}else{
		U = U0/((1-alpha)*L) * (L-x_red);
	}
	
	U *= mod;
	
	return U;
}

//----------------------------------------------------------
//					Definir la Fuerza
//----------------------------------------------------------
double MotorModel::force(double t, double tau_on, double tau_off){
	double tau = tau_on + tau_off;
	double F_cons;
	double mod; //modulador
	
	//Coordenadas reducidas
	double t_red = t - tau*floor(t/tau);
	double x_red = particle.getX() - L*floor(particle.getX()/L);
	
	//modulador
	if(0<=t_red && t_red<tau_on){
		mod = 1;
		
	}else{
		mod = 0;
	}
	
	double F;
	//fuerza
	if(0<x_red && x_red<alpha*L){
		F_cons = -U0/(alpha*L);
		
	}else if(alpha*L<x_red && x_red<L){
		F_cons = U0/((1-alpha)*L);
		
	}else{
		F_cons = 0;
	}
	
	F = mod*F_cons;
	
	return F;
}
