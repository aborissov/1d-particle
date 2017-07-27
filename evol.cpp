#include <iostream>
#include <chrono>
#include <random>
#include <math.h>
#include "constants.h"

using namespace std;

double calc_vel(double gamma, double beta){
        double v = sqrt(c*c - c*c/(gamma*gamma));
        //cout << "gamma = " << gamma << endl;
        //cout << "beta = " << beta << endl;
        //cout << "v = " << v << endl;
        return v/Vscl;
}

void move_particles(double *particles, int timestep, double *dw){
        double E,J,eta,eta_spitzer,nu,kappa,lambda_ei = 2.0e8/Lscl,Epar_extent = 1.0e3;
	double beta,v,u,uperp,gamma,dudt,betadot,gammadot,dbeta,dgamma,position;
	int random_index;
	random_index = 0;
	//random_index = timestep;

	eta_spitzer = 2.4e3/(pow(Temp,1.5))/etascl;
	J = 1.0e4/Escl;		// NON-DIMENSIONAL!!! Note: ensures electric field of 10 V/m when eta = 10^-3 (non-dimensional)
	eta = 1.0e-3;		// NON-DIMENSIONAL!!!
    	E = eta*J;		// NON-DIMENSIONAL!!!
	//cout << "Electric field: " << E*Escl << endl;


        for(int j = 0; j < nparticles; j++){
		position = particles[nfields*j];
		beta = particles[nfields*j+1];
		gamma = particles[nfields*j+2];
                v = calc_vel(gamma,beta);
                u = v*gamma*beta;
                uperp = v*gamma*sqrt(1.0-beta*beta);

		//printf("particle %d dw %f random_index %d\n",j,dw[random_index],random_index);

		if (abs(particles[nfields*j]*Lscl) < Epar_extent){
			//kappa = eta_spitzer/eta;
			kappa = 1.0e-5;
		}
		else {
			J = 0;
			E = 0;
			if (particles[nfields*j+3] == 0){
				particles[nfields*j+3] = timestep*dt*Tscl;
				//cout << "particle " << j << " " << timestep*dt*Tscl << endl;
			}
		}
		if (particles[nfields*j+3] == 0){
			if (abs(particles[nfields*j]*Lscl) < Epar_extent) nu = v/(lambda_ei*kappa);
			else nu = 0.0e6;
			//nu = 0.0;

                	dudt = q*E*Escl/m*Tscl/Vscl;
                	gammadot = u*Vscl/(c*c)*dudt*Vscl/(sqrt(1 + u*u*Vscl*Vscl/(c*c) + uperp*uperp*Vscl*Vscl/(c*c)));      // work in progress!!!
                	if (u == 0) betadot = 0;
                	else betadot = dudt/u*beta*(1.0-beta*beta);

                	dgamma = gammadot*dt;
                	dbeta = (betadot - beta*nu)*dt + sqrt((1.0 - beta*beta)*nu)*dw[random_index];
			random_index++;

                	particles[nfields*j+1] += dbeta;
                	if (particles[nfields*j+1] > 1.0){
						//cout << "beta = " << particles[nfields*j+1] << endl;
						particles[nfields*j+1] = -particles[nfields*j+1] + floor(particles[nfields*j+1]) + 1.0;
					}
                	else if (particles[nfields*j+1] < -1.0){
						//cout << "beta = " << particles[nfields*j+1] << endl;
						particles[nfields*j+1] = -particles[nfields*j+1] + ceil(particles[nfields*j+1]) - 1.0;
					}
                	particles[nfields*j+2] += dgamma;
                	if (particles[nfields*j+2] < 1) particles[nfields*j+2] -= 2.0*dgamma;

                	v = calc_vel(particles[nfields*j + 2], particles[nfields*j + 1]);
                	particles[nfields*j] += particles[nfields*j+1]*v*dt;
			
			energy_kev[j] = (particles[nfields*j+2]-1)*511.0;
			potential[j] = -eta*J*Escl*particles[nfields*j]*Lscl/1.0e3;
			

			if (fabs(energy_kev[j] - potential[j] - energy_kev_0) > 1.0){
				particles[nfields*j] = Epar_extent/Lscl;
				printf("particle %d deviated from energy conservation at time %f", j, timestep*dt*Tscl);
				printf(" kinetic %f, potential %f, initial %f, difference %f\n", energy_kev[j], potential[j], energy_kev_0,energy_kev[j] - potential[j] - energy_kev_0);
			}
			//cout << "particle " << j << " total energy " << energy_kev[j] - potential[j] - energy_kev_0  
			//	<< " kinetic, potential, initial " << energy_kev[j] << " " << potential[j] << " "  << energy_kev_0 << endl;
		}
        }
}

