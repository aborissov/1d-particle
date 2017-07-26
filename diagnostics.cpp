#include <iostream>
#include <fstream>
#include "constants.h"
#include "evol.h"

using namespace std;

void write_particles(double *particles, bool newflag, char *filepath){
        ofstream outFile;
        double v;
        v = calc_vel(particles[2], particles[1]);
        double pos,theta,energy,t_final;

        if (newflag) outFile.open(filepath, ofstream::binary);
        else outFile.open(filepath, ofstream::binary | ofstream::app);
	outFile.write((char*) &(nparticles), sizeof(int));
	outFile.write((char*) &(output_fields), sizeof(int));
	for (int j = 0; j < nparticles; j++){
		pos = particles[nfields*j]*Lscl;
        	outFile.write((char*) &(pos), sizeof(double));
	}
	for (int j = 0; j < nparticles; j++){
		energy = (particles[nfields*j+2]-1)*5.11e5;
        	outFile.write((char*) &(energy), sizeof(double));
		//cout << "energy = " << energy << endl;
	}
	for (int j = 0; j < nparticles; j++){
		theta = particles[nfields*j+1];
        	outFile.write((char*) &(theta), sizeof(double));
	}
	for (int j = 0; j < nparticles; j++){
		t_final = particles[nfields*j+3]*Tscl;
        	outFile.write((char*) &(t_final), sizeof(double));
	}
        outFile.close();
}

void write_particle(double *particles, bool newflag){
	ofstream outFile;
	double v = calc_vel(particles[2], particles[1]);
	double pos,theta,energy,t_final;

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
		outFile.write((char *) &pos, sizeof(double));
		outFile.write((char *) &theta, sizeof(double));
		outFile.write((char *) &energy, sizeof(double));
		outFile.write((char *) &t_final, sizeof(double));
	}

}
