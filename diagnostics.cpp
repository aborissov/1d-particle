#include <iostream>
#include <fstream>
#include "constants.h"
#include "evol.h"

using namespace std;

void write_particles(float *particles, bool newflag, char *filepath){
        ofstream outFile;
        float v;
        v = calc_vel(particles[2], particles[1]);
        float pos,theta,energy;

        if (newflag) outFile.open(filepath, ofstream::binary);
        else outFile.open(filepath, ofstream::binary | ofstream::app);
	outFile.write((char*) &(nparticles), sizeof(int));
	outFile.write((char*) &(output_fields), sizeof(int));
	for (int j = 0; j < nparticles; j++){
		pos = particles[nfields*j]*Lscl;
        	outFile.write((char*) &(pos), sizeof(float));
	}
	for (int j = 0; j < nparticles; j++){
		energy = (particles[nfields*j+2]-1)*5.11e5;
        	outFile.write((char*) &(energy), sizeof(float));
		//cout << "energy = " << energy << endl;
	}
	for (int j = 0; j < nparticles; j++){
		theta = particles[nfields*j+1];
        	outFile.write((char*) &(theta), sizeof(float));
	}
        //outFile.write((char*) &v, sizeof(float));
        outFile.close();
        //cout << particles[0] << endl;
}

void write_particle(float *particles, bool newflag){
	ofstream outFile;
	float v = calc_vel(particles[2], particles[1]);
	float pos,theta,energy;

	if (newflag){
		outFile.open("trajectories.dat", ofstream::binary);
		outFile.write((char *) &(nparticles), sizeof(int));
		outFile.write((char *) &(output_fields), sizeof(int));
		outFile.write((char *) &(nwrite), sizeof(int));
	}
	else outFile.open("trajectories.dat", ofstream::binary | ofstream::app);
	for (int j = 0; j < nwrite_particles; j++){
		pos = particles[nfields*j]*Lscl;
		theta = particles[nfields*j+1];
		energy = (particles[nfields*j+2] - 1)*5.11e5;
		outFile.write((char *) &pos, sizeof(float));
		outFile.write((char *) &theta, sizeof(float));
		outFile.write((char *) &energy, sizeof(float));
	}

}