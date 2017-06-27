#include <iostream>
#include <chrono>
#include <random>
#include <math.h>
#include "constants.h"

using namespace std;

float calc_vel(float gamma, float beta){
        float v = sqrt(c*c - c*c/(gamma*gamma));
        //cout << "gamma = " << gamma << endl;
        //cout << "beta = " << beta << endl;
        //cout << "v = " << v << endl;
        return v/Vscl;
}

void move_particles(float *particles, int timestep){
        float q,E,m,J,eta,eta_spitzer,Temp,nu,kappa,lambda_ei = 2.0e8/Lscl,Epar_extent = 2.0e4;
        q = -1.6e-19;
        m = 9.11e-31;
	Temp = 1.0e7;
	eta_spitzer = 2.4e3/(pow(Temp,1.5))/etascl;
	J = 10.0;	// NON-DIMENSIONAL!!!
	eta = 1.0e-5;	// NON-DIMENSIONAL!!!
        E = eta*J;	// NON-DIMENSIONAL!!!


        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator (seed);
        for(int j = 0; j < nparticles; j++){
		if (abs(particles[nfields*j]*Lscl) < Epar_extent){
			kappa = eta_spitzer/eta;
		}
		else {
			J = 0;
			E = 0;
			particles[nfields*j+3] = timestep*dt*Tscl;
		}
		if (particles[nfields*j+3] == 0){
                //unsigned seed = chrono::system_clock::now().time_since_epoch().count();
                //default_random_engine generator (seed*(j+1));
                normal_distribution<float> distribution(0.0,1.0);
                float dw = distribution(generator)*sqrt(dt);

                float beta = particles[nfields*j+1];
                float v = calc_vel(particles[nfields*j+2], particles[nfields*j+1]);
				if (abs(particles[nfields*j]*Lscl) < Epar_extent) nu = v/(lambda_ei*kappa);
				//if (-1 > 0) nu = v/(lambda_ei*kappa);
				else nu = 0.0e6;
                float u = v*particles[nfields*j+2]*beta;
                float uperp = v*particles[nfields*j+2]*sqrt(1.0-beta*beta);
                float dudt = q*E*Escl/m*Tscl/Vscl;
                float gammadot = u*Vscl/(c*c)*dudt*Vscl/(sqrt(1 + u*u*Vscl*Vscl/(c*c) + uperp*uperp*Vscl*Vscl/(c*c)));      // work in progress!!!
                float betadot;
                if (u == 0) betadot = 0;
                else betadot = dudt/u*beta*(1.0-beta*beta);

                float dgamma = gammadot*dt;
                float dbeta = (betadot - beta*nu)*dt + sqrt((1.0 - beta*beta)*nu)*dw;

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
                //particles[nfields*j] = dw/sqrt(dt)/Lscl;
                //particles[nfields*j+2] = sqrt(1.0 + u*u/(c*c) + 2.0*mu*B/(m*c*c));

                //cout << "z = " << particles[0]*Lscl << endl;
                //cout << "beta = " << beta << endl;
                //cout << "gamma = " << particles[2] << endl;
                //cout << "v = " << v*Vscl << endl;
                //cout << "u = " << u << endl;
                //cout << "uperp = " << uperp << endl;
                //cout << "dudt = " << dudt << endl;
                //cout << "gammadot = " << gammadot << endl;
                //cout << "betadot = " << betadot << endl;
		}
        }
}

