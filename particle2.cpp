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
	double particles[nparticles*nfields]; // array of particles: 1d position, cosine of pitch angle, lorentz factor, exit time
	bool newflag = 1,newflag_trajectories = 1;
	double maxp = 0,minp = 0,meanp = 0;
	double *dw;

	energy_kev = (double *)malloc(nparticles*sizeof(double));
	potential = (double *)malloc(nparticles*sizeof(double));
	dw = (double *)malloc(nparticles*nt*sizeof(double));
	initialise(particles);

	// initialise random number generator
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator (seed);
	normal_distribution<double> distribution(0.0,1.0);

	// main time loop
	for (int j = 0; j < nt; j++){
		//if (j%1 == 0) write_particle(particles,newflag_trajectories);
		//newflag_trajectories = 0;

		// generate random numbers
		for (int j = 0; j < nparticles; j++) dw[j] = distribution(generator)*sqrt(dt);

		move_particles(particles,j,dw);
		if (isnan(particles[0])){
			cout << "position is nan. stopping" << endl;
			return 0;
		}
	}

	free(energy_kev);
	free(potential);
	free(dw);

	write_particles(particles,newflag,argv[1]);
	cout << "size of particles array " << nparticles << endl;
	for (int j = 0; j < nparticles; j++) printf("particle %d final energy %f position %f\n",j,(particles[nfields*j+2]-1)*511.0,particles[nfields*j]*Lscl);

	return 0;
}
