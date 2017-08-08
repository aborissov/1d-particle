#include <iostream>
#include <chrono>
#include <random>
#include <math.h>
#include <omp.h>
#include "constants.h"

using namespace std;

real calc_vel(real gamma, real beta){
        real v = sqrt(c*c - c*c/(gamma*gamma));
        //cout << "gamma = " << gamma << endl;
        //cout << "beta = " << beta << endl;
        //cout << "v = " << v << endl;
        return v/Vscl;
}

void move_particles(real *particles, int timestep, real *dw){
        real E,J,eta,eta_spitzer,nu,kappa,lambda_ei = 2.0e8/Lscl,Epar_extent = 1.0e3;
		real beta,v,u,uperp,gamma,dudt,betadot,gammadot,dbeta,dgamma,position;
		int random_index;
		random_index = 0;
		//random_index = timestep;

		int np,tid;
		omp_set_dynamic(0);
		omp_set_num_threads(nparticles);
		#pragma omp parallel private(np,position,beta,gamma,v,u,uperp,random_index,tid,dudt,gammadot,betadot,dgamma,dbeta,J,E)
		{
		#ifdef _OPENMP
		tid = omp_get_thread_num();
		random_index = tid;
		np = tid;
		#endif

        //for(np = 0; np < nparticles; np++){
		position = particles[nfields*np];
		beta = particles[nfields*np+1];
		gamma = particles[nfields*np+2];
        v = calc_vel(gamma,beta);
        u = v*gamma*beta;
        uperp = v*gamma*sqrt(1.0-beta*beta);

		//#pragma omp critical
		//{
		//printf("thread %d particle %d dw %f random_index %d timestep %d\n",tid,np,dw[random_index],random_index,timestep);
		//}

		if (abs(position*Lscl) > Epar_extent){
			J = 0;
			E = 0;
			if (particles[nfields*np+3] == 0){
				particles[nfields*np+3] = timestep*dt;
				//cout << "particle " << np << " " << timestep*dt*Tscl << endl;
			}
		}
		else{
		if (particles[nfields*np+3] == 0){
			eta_spitzer = 2.4e3/(pow(Temp,1.5))/etascl;
			J = 1.0e4/Escl;		// NON-DIMENSIONAL!!! Note: ensures electric field of 10 V/m when eta = 10^-3 (non-dimensional)
			eta = 1.0e-3;		// NON-DIMENSIONAL!!!
    		E = eta*J;		// NON-DIMENSIONAL!!!
			//kappa = eta_spitzer/eta;
			kappa = 1.0e-5;

			if (abs(position*Lscl) < Epar_extent) nu = v/(lambda_ei*kappa);
			else nu = 0.0;
			//nu = 0.0;

            dudt = q*E*Escl/m*Tscl/Vscl;
            gammadot = u*Vscl/(c*c)*dudt*Vscl/(sqrt(1 + u*u*Vscl*Vscl/(c*c) + uperp*uperp*Vscl*Vscl/(c*c)));      // work in progress!!!
            if (u == 0) betadot = 0;
            else betadot = dudt/u*beta*(1.0-beta*beta);

            dgamma = gammadot*dt;
            dbeta = (betadot - beta*nu)*dt + sqrt((1.0 - beta*beta)*nu)*dw[random_index];
			random_index++;
			#ifdef _OPENMP
			random_index = tid;
			#endif

            beta += dbeta;
            if (beta > 1.0){
				//cout << "beta = " << particles[nfields*np+1] << endl;
				beta = -beta + floor(beta) + 1.0;
			}
            else if (beta < -1.0){
				//cout << "beta = " << particles[nfields*np+1] << endl;
				beta = -beta + ceil(beta) - 1.0;
			}
            gamma += dgamma;
            if (gamma < 1) gamma -= 2.0*dgamma;

            v = calc_vel(gamma, beta);
            position += beta*v*dt;
			
			energy_kev[np] = (gamma-1.0)*511.0;
			potential[np] = -eta*J*Escl*position*Lscl/1.0e3;
			

			particles[nfields*np] = position;
			particles[nfields*np+1] = beta;
			particles[nfields*np+2] = gamma;
			
			//if (fabs(energy_kev[np] - potential[np] - energy_kev_0) < 1.0 && np == 1){
			//	printf("particle %d deviated from energy conservation at time %f", np, timestep*dt*Tscl);
			//	printf(" kinetic %f, potential %f, initial %f, difference %f\n", energy_kev[np], potential[np], energy_kev_0,energy_kev[np] - potential[np] - energy_kev_0);
			//	printf(" eta, J, Escl, position, epar_extent, %f %f %f %f %f\n", eta,J,Escl,position*Lscl,Epar_extent);
			//}
			if (fabs(energy_kev[np] - potential[np] - energy_kev_0) > 1.0){
				particles[nfields*np] = Epar_extent/Lscl;
				printf("particle %d deviated from energy conservation at time %f", np, timestep*dt*Tscl);
				printf(" kinetic %f, potential %f, initial %f, difference %f\n", energy_kev[np], potential[np], energy_kev_0,energy_kev[np] - potential[np] - energy_kev_0);
				printf(" eta, J, Escl, position, epar_extent, %f %f %f %f %f\n", eta,J,Escl,position*Lscl,Epar_extent);
			}
			//cout << "particle " << np << " total energy " << energy_kev[np] - potential[np] - energy_kev_0  
			//	<< " kinetic, potential, initial " << energy_kev[np] << " " << potential[np] << " "  << energy_kev_0 << endl;
		}
        }
		}
}

