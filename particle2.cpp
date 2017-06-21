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

int main(){
	float particles[nparticles*nfields]; // array of particles: 1d position, cosine of pitch angle, lorentz factor
	bool newflag = 1;

	initialise(particles,nparticles);
	for (int j = 0; j < nt; j++){
		//if (j%1 == 0) write_particles(particles,nparticles,newflag);
		//newflag = 0;
		move_particles(particles,nparticles);
		//cout << particles[0] << endl;
		//cout << "time = " << j*dt*Tscl << endl;
		if (isnan(particles[0])){
			cout << "position is nan. stopping" << endl;
			return 0;
		}
	}

	write_particles(particles,newflag);
	cout << "size of particles array " << nparticles << endl;
	return 0;
}
