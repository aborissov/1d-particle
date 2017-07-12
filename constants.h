#ifndef CONSTANTS_H
#define CONSTANTS_H

const float c = 3.0e8;
const float mu0_si = 8.85e-12;

const float Bscl = 0.01;
const float Lscl = 1.0e6;
const float Tscl = 1.0e2;
const float Vscl = Lscl/Tscl;
const float Escl = Vscl*Bscl;
const float etascl = Vscl*Lscl*mu0_si;

const float Tfinal = 0.1/Tscl;
const float dt = 1.0e-7/Tscl;
const int nt = Tfinal/dt;
const int nfields = 4;
const int nparticles = 10000;
const int nwrite = 500;
const int nwrite_particles = 10;
const int output_fields = nfields;

#endif
