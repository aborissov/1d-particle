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
        double q,E,m,J,eta,eta_spitzer,Temp,nu,kappa,lambda_ei = 2.0e8/Lscl,Epar_extent = 1.0e3;
	int random_index = 0;

        q = -1.6e-19;
        m = 9.11e-31;
		Temp = 1.0e7;
		eta_spitzer = 2.4e3/(pow(Temp,1.5))/etascl;
		J = 10.0;		// NON-DIMENSIONAL!!!
		eta = 1.0e-2;		// NON-DIMENSIONAL!!!
    		E = eta*J;		// NON-DIMENSIONAL!!!
		//cout << "Electric field: " << E*Escl << endl;


        //unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        //default_random_engine generator (seed);
        for(int j = 0; j < nparticles; j++){
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
                //unsigned seed = chrono::system_clock::now().time_since_epoch().count();
                //default_random_engine generator (seed*(j+1));
                //normal_distribution<double> distribution(0.0,1.0);
                //double dw = distribution(generator)*sqrt(dt);

                double beta = particles[nfields*j+1];
                double v = calc_vel(particles[nfields*j+2], particles[nfields*j+1]);

				if (abs(particles[nfields*j]*Lscl) < Epar_extent) nu = v/(lambda_ei*kappa);
				else nu = 0.0e6;
				//nu = 0.0;

                double u = v*particles[nfields*j+2]*beta;
                double uperp = v*particles[nfields*j+2]*sqrt(1.0-beta*beta);
                double dudt = q*E*Escl/m*Tscl/Vscl;
                double gammadot = u*Vscl/(c*c)*dudt*Vscl/(sqrt(1 + u*u*Vscl*Vscl/(c*c) + uperp*uperp*Vscl*Vscl/(c*c)));      // work in progress!!!
                double betadot;
                if (u == 0) betadot = 0;
                else betadot = dudt/u*beta*(1.0-beta*beta);

                double dgamma = gammadot*dt;
                double dbeta = (betadot - beta*nu)*dt + sqrt((1.0 - beta*beta)*nu)*dw[random_index];
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
		
		//energy_kev[j] = (particles[nfields*j+2]-1)*511.0;
		//potential[j] = -10.0*particles[nfields*j]*Lscl/1.0e3;
		//cout << "particle " << j << " total energy " << energy_kev[j] - potential[j] - energy_kev_0  
		//	<< " kinetic, potential, initial " << energy_kev[j] << " " << potential[j] << " "  << energy_kev_0 << endl;
		}
        }
}

