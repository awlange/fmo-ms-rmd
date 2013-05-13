/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "dynamics.h"
#include "atom.h"
#include "run.h"
#include "matrix.h"
#include "input.h"

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Run Velocity Verlet molecular dynamics 
  The big outer loop over MD steps
-----------------------------------------------------------------*/
void Dynamics::runMolecularDynamics()
{
  // Before calling this function, one must have already init() 

  // ** Reporting for step zero ** //
  stepReport(0);

  // *** The loop *** //
  for (int istep=0; istep<nTimeSteps; ++istep) {

    // Increment global time step index, which may have been read in
    iCurrentStep++;
    
    // ** Update the position ** //
    // r(t + dt) = r(t) + v(t)*dt + 0.5*a(t)*dt*dt
    // Note: a(t) = -(1/m)*dV( r(t) )
    stepPosition();

    // ** Update the half-step velocity ** //
    // v(t + 0.5*dt) = v(t) + 0.5*a(t)*dt
    halfStepVelocity(0); 

    // ** Update the acceleration (i.e. force) ** //
    // Need to recompute the force with the new coordinates
    // a(t + dt) = -(1/m)*dV( r(t + dt) )
    fmr->run->calculate_force();
    
    // ** Update the full-step velocity ** //
    // v(t + dt) = v(t + 0.5*dt) + 0.5*a(t + dt)*dt
    halfStepVelocity(1); 

    // ** Remove translation and rotation from velocity ** //
    if (iThermostat != THERMOSTAT_OFF) {
      removeTransRotFromVelocity();
    }

    // ** Handle thermostating ** //
    applyThermostat();

    // ** Reporting for step ** //
    stepReport(istep+1);
  }
}

/*-----------------------------------------------------------------
  Update the position
-----------------------------------------------------------------*/
void Dynamics::stepPosition()
{
  // Master does calculation, broadcasts to rest
  int natoms = fmr->atom->natoms;
  double *coord = fmr->atom->coord;
  if (fmr->master_rank) {
    double *veloc = fmr->atom->veloc;
    double *force = fmr->atom->force;
    double *mass  = fmr->atom->mass;
    double hdt2   = 0.5*dt*dt;
    double b2a    = fmr->math->b2a;

    // Move the atom positions forward, making sure to convert to Angstrom
    for (int i=0; i<natoms; ++i) {
      coord[3*i]   += (veloc[3*i  ]*dt - force[3*i  ]*hdt2/mass[i]) * b2a;
      coord[3*i+1] += (veloc[3*i+1]*dt - force[3*i+1]*hdt2/mass[i]) * b2a;
      coord[3*i+2] += (veloc[3*i+2]*dt - force[3*i+2]*hdt2/mass[i]) * b2a;
    }
  }
  MPI_Bcast(coord, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
}

/*-----------------------------------------------------------------
  Update the velocity half a step 
-----------------------------------------------------------------*/
void Dynamics::halfStepVelocity(int remove)
{
  // Master does calculation, broadcasts to rest
  int natoms = fmr->atom->natoms;
  double *veloc = fmr->atom->veloc;
  if (fmr->master_rank) {
    double *force = fmr->atom->force;
    double *mass  = fmr->atom->mass;
    double hdt    = 0.5*dt;
  
    for (int i=0; i<natoms; ++i) {
      veloc[3*i  ] -= force[3*i  ]*hdt/mass[i];
      veloc[3*i+1] -= force[3*i+1]*hdt/mass[i];
      veloc[3*i+2] -= force[3*i+2]*hdt/mass[i];
    }
  }
  MPI_Bcast(veloc, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
}


/*-----------------------------------------------------------------
  Report information about MD for this step 
-----------------------------------------------------------------*/
void Dynamics::stepReport(int istep)
{
  EPotential = fmr->matrix->GSEnergy;
  EKinetic   = computeKE();
  ETotal     = EPotential + EKinetic;
  currentTemperature = computeTemperature();

  if (fmr->master_rank) {
    double au2kcal = fmr->math->au2kcal;
    printf("---------------------------------------------------------------------\n");
    printf("STEP: %d\n", iCurrentStep);
    printf("EPotential  = %16.8f hartree = %16.8f kcal/mol\n", EPotential, EPotential*au2kcal);
    printf("EKinetic    = %16.8f hartree = %16.8f kcal/mol\n", EKinetic, EKinetic*au2kcal);
    printf("ETotal      = %16.8f hartree = %16.8f kcal/mol\n", ETotal, ETotal*au2kcal);
    printf("Temperature = %10.6f Kelvin\n", currentTemperature);
    printf("---------------------------------------------------------------------\n");
  }

  // ** Write coordinates to file ** //
  writeTrajCoords(istep);

  // ** Write restart file ** //
  fmr->input->write_restart_file();
}
