/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "dynamics.h"
#include "atom.h"
#include "math.h"

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Set atom velocities based on vOption 
-----------------------------------------------------------------*/
void Dynamics::setVelocities()
{
  int natoms    = fmr->atom->natoms;
  double *veloc = fmr->atom->veloc;

  if (vOption == DYNAMICS_ZERO) {
    // ** Set all velocities to zero ** //
    for (int i=0; i<natoms; ++i) {
      veloc[3*i]   = 0.0; 
      veloc[3*i+1] = 0.0; 
      veloc[3*i+2] = 0.0; 
    } 
  }
  else if (vOption == DYNAMICS_READ) {
    // Do nothing. These should have been set in read_restart_file. 
  }
  else if (vOption == DYNAMICS_RAND) {
    // ** Random velocity directions with temperature rescaling ** //
  
    if (fmr->master_rank) {

      // Seed if needed
      fmr->math->rng_seed(vseed);

      // Choose a random velocity 
      for (int i=0; i<natoms; ++i) { 
        veloc[3*i]   = fmr->math->rngInRange(-1.0, 1.0);
        veloc[3*i+1] = fmr->math->rngInRange(-1.0, 1.0);
        veloc[3*i+2] = fmr->math->rngInRange(-1.0, 1.0);
      }
    }

    // Clean out translation and rotations from the above
    removeTransRotFromVelocity();
  
    if (fmr->master_rank) {
      // Scale velocities in accord with target temperature
      double T = computeTemperature(); 
      double lambda = sqrt(targetTemperature / T);
      for (int i=0; i<natoms; ++i) {
        veloc[3*i]   *= lambda;
        veloc[3*i+1] *= lambda;
        veloc[3*i+2] *= lambda;
      }
    }
    // Send out the scaled velocities to others
    MPI_Bcast(veloc, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world); 
  }

  // ** Printing ** //
  if (fmr->master_rank) {
    printf("Initial velocities (in a.u.):\n");
    for (int i=0; i<natoms; ++i) {
      printf("%3d %14.8f %14.8f %14.8f\n", i, veloc[3*i], veloc[3*i+1], veloc[3*i+2]);
    }
  }
}


