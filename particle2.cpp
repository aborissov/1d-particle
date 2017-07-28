#define MAINFILE
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <chrono>
#include "constants.h"
#include "initial_conditions.h"
#include "evol.h"
#include "diagnostics.h"

using namespace std;

int main(int argc, char *argv[]){
	real particles[nparticles*nfields]; // array of particles: 1d position, cosine of pitch angle, lorentz factor, exit time
	bool newflag = 1,newflag_trajectories = 1;
	real *dw;

	energy_kev = (real *)malloc(nparticles*sizeof(real));
	potential = (real *)malloc(nparticles*sizeof(real));
	//dw = (real *)malloc(nparticles*sizeof(real));
	dw = (real *)malloc(nt*sizeof(real));
	initialise(particles);

	// initialise random number generator
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	//unsigned seed1 = 1234;
	//unsigned seed2 = 29485;
	//default_random_engine generator1 (seed1);
	//default_random_engine generator2 (seed2);
	default_random_engine generator (seed);
	normal_distribution<real> distribution(0.0,1.0);
	//for (int j = 0; j < nt; j++) dw[j] = distribution(generator)*sqrt(dt);

	// main time loop
	for (int j = 0; j < nt; j++){
		//if (j%1 == 0) write_particle(particles,newflag_trajectories);
		//newflag_trajectories = 0;

		// generate random numbers
		for (int j = 0; j < nparticles; j++) dw[j] = distribution(generator)*sqrt(dt);
		//dw[1] = distribution(generator1)*sqrt(dt);
		//dw[0] = distribution(generator2)*sqrt(dt);
		//dw[0] = dw[1];

		move_particles(particles,j,dw);
		if (isnan(particles[0])){
			cout << "position is nan. stopping" << endl;
			return 0;
		}
		//if (j % (nt/100) == 0) printf("timestep %d of %d\n",j,nt);
	}

	free(energy_kev);
	free(potential);
	free(dw);

	write_particles(particles,newflag,argv[1]);
	cout << "size of particles array " << nparticles << endl;
	//for (int j = 0; j < nparticles; j++) printf("particle %d final energy %f position %f\n",j,(particles[nfields*j+2]-1)*511.0,particles[nfields*j]*Lscl);

	return 0;
}
