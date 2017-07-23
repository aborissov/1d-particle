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
	float particles[nparticles*nfields]; // array of particles: 1d position, cosine of pitch angle, lorentz factor
	bool newflag = 1,newflag_trajectories = 1;
	float maxp = 0,minp = 0,meanp = 0;

	initialise(particles);
	for (int j = 0; j < nt; j++){
		if (j%1 == 0) write_particle(particles,newflag_trajectories);
		newflag_trajectories = 0;
		move_particles(particles,j);
		if (isnan(particles[0])){
			cout << "position is nan. stopping" << endl;
			return 0;
		}
	}

	write_particles(particles,newflag,argv[1]);
	cout << "size of particles array " << nparticles << endl;

	return 0;
}
