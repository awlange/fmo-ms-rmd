/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "dynamics.h"
#include "atom.h"
#include "run.h"
#include "math.h"

using namespace FMR_NS;

/*-----------------------------------------------------------------
  General thermostat function for running NVT
-----------------------------------------------------------------*/
void Dynamics::applyThermostat()
{
  // ** Call the selected thermostat option ** //
  if (iThermostat == THERMOSTAT_OFF) {
    // No thermostating. Running NVE.
    return;
  } else if (iThermostat == THERMOSTAT_LANGEVIN) {
    thermostatLangevin();
  } else if (iThermostat == THERMOSTAT_BERENDSEN) {
    thermostatBerendsen();
  }
}

/*-----------------------------------------------------------------
  Langevin thermostat 
-----------------------------------------------------------------*/
void Dynamics::thermostatLangevin()
{
  double *force = fmr->atom->force;
  double *veloc = fmr->atom->veloc;
  double *mass  = fmr->atom->mass;
  int natoms    = fmr->atom->natoms;

  if (fmr->master_rank) {
    // Check seeding. We may not have seeded yet if not using random velocities to start.
    fmr->math->rng_seed(vseed);

    double gamma = 1.0/tau;
    double factor = 2.0*gamma*(fmr->math->kB_au)*targetTemperature / dt;
    for (int i=0; i < natoms; ++i) {
      const double fm  = sqrt(factor * mass[i]);
      const double rfx = fm * fmr->math->rngInRange(-1.0, 1.0); 
      const double rfy = fm * fmr->math->rngInRange(-1.0, 1.0);
      const double rfz = fm * fmr->math->rngInRange(-1.0, 1.0);
      const double gm  = gamma * mass[i];
      force[3*i  ] += -gm*veloc[3*i  ] + rfx;
      force[3*i+1] += -gm*veloc[3*i+1] + rfy;
      force[3*i+2] += -gm*veloc[3*i+2] + rfz;
    }
  }

  // ** Master broadcasts adjusted forces here ** //
  MPI_Bcast(force, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
}

/*-----------------------------------------------------------------
  Berendsen thermostat 
-----------------------------------------------------------------*/
void Dynamics::thermostatBerendsen()
{
  if (fmr->master_rank) {
    currentTemperature = computeTemperature();
    double lambda = sqrt(1.0 + (dt/tau) * (targetTemperature/currentTemperature - 1.0) );
    for (int i=0; i < fmr->atom->natoms; ++i) {
      fmr->atom->veloc[3*i]   *= lambda;
      fmr->atom->veloc[3*i+1] *= lambda;
      fmr->atom->veloc[3*i+2] *= lambda;
    }
  }
  // ** Master broadcasts adjusted velocities here ** //
  MPI_Bcast(fmr->atom->veloc, 3*fmr->atom->natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
}
