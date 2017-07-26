#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifdef MAINFILE
	#define EXTERN
#else
	#define EXTERN extern
#endif

const double c = 3.0e8;
const double mu0_si = 8.85e-12;

const double Bscl = 0.01;
const double Lscl = 1.0e6;
const double Tscl = 1.0e2;
const double Vscl = Lscl/Tscl;
const double Escl = Vscl*Bscl;
const double etascl = Vscl*Lscl*mu0_si;

const double Tfinal = 0.1/Tscl;
const double dt = 1.0e-8/Tscl;
const int nt = Tfinal/dt;
const int nfields = 4;
const int nparticles = 5;
const int nwrite = 500;
const int nwrite_particles = 1;
const int output_fields = nfields;

EXTERN double *energy_kev, *potential;
const double energy_kev_0 = 0.5*511.0;

#endif
