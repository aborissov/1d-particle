#include <iostream>
#include "constants.h"

using namespace std;

void initialise(float *particles,int nparticles){
        for(int j = 0; j < nparticles*nfields; j++){
                if (j%nfields == 0) particles[j] = 0;
                else if (j%nfields == 1) particles[j] = -0.1;
                else particles[j] = 1.5;
        }
}
