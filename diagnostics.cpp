#include <iostream>
#include <fstream>
#include "constants.h"
#include "evol.h"

using namespace std;

void write_particles(real *particles, bool newflag, char *filepath){
        ofstream outFile;
        real v;
        v = calc_vel(particles[2], particles[1]);
        real pos,theta,energy,t_final;

        if (newflag) outFile.open(filepath, ofstream::binary);
        else outFile.open(filepath, ofstream::binary | ofstream::app);
	outFile.write((char*) &(nparticles), sizeof(int));
	outFile.write((char*) &(output_fields), sizeof(int));
	for (int j = 0; j < nparticles; j++){
		pos = particles[nfields*j]*Lscl;
        	outFile.write((char*) &(pos), sizeof(real));
	}
	for (int j = 0; j < nparticles; j++){
		energy = (particles[nfields*j+2]-1)*5.11e5;
        	outFile.write((char*) &(energy), sizeof(real));
		//cout << "energy = " << energy << endl;
	}
	for (int j = 0; j < nparticles; j++){
		theta = particles[nfields*j+1];
        	outFile.write((char*) &(theta), sizeof(real));
	}
	for (int j = 0; j < nparticles; j++){
		t_final = particles[nfields*j+3]*Tscl;
        	outFile.write((char*) &(t_final), sizeof(real));
	}
        outFile.close();
}

void write_particle(real *particles, bool newflag){
	ofstream outFile;
	real v = calc_vel(particles[2], particles[1]);
	real pos,theta,energy,t_final;

	if (newflag){
		outFile.open("trajectories.dat", ofstream::binary);
		outFile.write((char *) &(nparticles), sizeof(int));
		outFile.write((char *) &(output_fields), sizeof(int));
		outFile.write((char *) &(nwrite_particles), sizeof(int));
	}
	else outFile.open("trajectories.dat", ofstream::binary | ofstream::app);
	for (int j = 0; j < nwrite_particles; j++){
		pos = particles[nfields*j]*Lscl;
		theta = particles[nfields*j+1];
		energy = (particles[nfields*j+2] - 1)*5.11e5;
		t_final = particles[nfields*j+3]*Tscl;
		outFile.write((char *) &pos, sizeof(real));
		outFile.write((char *) &theta, sizeof(real));
		outFile.write((char *) &energy, sizeof(real));
		outFile.write((char *) &t_final, sizeof(real));
	}

}
