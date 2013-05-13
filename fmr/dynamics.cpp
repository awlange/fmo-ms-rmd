/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "dynamics.h"
#include "state.h"

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/

Dynamics::Dynamics(FMR *fmr) : Pointers(fmr)
{
  // Basic initializations
  nTimeSteps         = 0;
  iCurrentStep       = 0;
  iThermostat        = THERMOSTAT_OFF;
  dt                 = 0.0;
  vOption            = 0;
  targetTemperature  = 0.0;
  currentTemperature = 0.0;
  EPotential         = 0.0;
  EKinetic           = 0.0;
  ETotal             = 0.0;
  vseed              = -1; // deafult to using timestamp, read in otherwise if desired
  sprintf(trajFile, "%s", "fmr_traj.xyz"); // default file name 
  tau                = 0.5; // default, femtoseconds
}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/

Dynamics::~Dynamics()
{

}


/*-----------------------------------------------------------------
  Initializer 
-----------------------------------------------------------------*/
void Dynamics::init()
{
  // Note: input.cpp read in various MD run parameters already before this point
  //       Also, this is to be called after a RUN_FORCE calculation is completed

  if (fmr->master_rank) {
    printf("Preparing to run molecular dynamics for %d timesteps.\n", nTimeSteps);
  }

  // ** Initialize the atom velocities ** // 
  // Also, remove translation and rotation from velocity and gradient
  setVelocities();
  removeTransRotFromVelocity();
  //removeTransRotFromGradient();

  // ** Initialize any stats we want to keep ** //

  // ** Other miscellaneous business ** //
  // Make subsequent SCF calls read from disk, helps to speed a little
  //fmr->state->flag_read_MOs = 1;
}


/*-----------------------------------------------------------------
  Compute the current kinetic energy 
-----------------------------------------------------------------*/
double Dynamics::computeKE() 
{
  int natoms    = fmr->atom->natoms;
  double *veloc = fmr->atom->veloc;
  double *mass  = fmr->atom->mass;
  double KE = 0.0;
  for (int i=0; i<natoms; ++i) {
    const double tmp = veloc[3*i]   * veloc[3*i] +
                       veloc[3*i+1] * veloc[3*i+1] +
                       veloc[3*i+2] * veloc[3*i+2];
    KE += tmp * mass[i];
  }
  return KE *= 0.5;
}

/*-----------------------------------------------------------------
  Compute the current temperature 
-----------------------------------------------------------------*/
double Dynamics::computeTemperature() 
{
  return computeKE() / ( (3.0/2.0)*fmr->math->kB_au*((double)fmr->atom->natoms) );
}


/*-----------------------------------------------------------------
  Write the current positions and velocities to the trajectory file 
  mode 0 = write new file
  any other mode = append to current file
  same as write restart file, but here we append
-----------------------------------------------------------------*/
void Dynamics::writeTrajCoords(int mode) 
{
  if (fmr->master_rank) {
    int natoms    = fmr->atom->natoms;
    double *coord = fmr->atom->coord;
    double *veloc = fmr->atom->veloc;
    FILE *fs;

    // Check mode
    if (mode == 0) fs = fopen(trajFile, "w");
    else           fs = fopen(trajFile, "a");
    if (fs == NULL) {
      char tmpstr[256];
      sprintf(tmpstr, "Failure to write to trajectory file %s", trajFile);
      fmr->error(FLERR, tmpstr);
    }

    fprintf(fs, "%d %d %d\n", natoms, fmr->atom->nfragments, fmr->atom->ireactive);
    fprintf(fs, "%d\n", fmr->dynamics->iCurrentStep);
    for (int i=0; i<natoms; ++i) {
      fprintf(fs, "%c %16.12f %16.12f %16.12f %3d %16.12f %16.12f %16.12f\n",
                  fmr->atom->symbol[i],
                  coord[3*i], coord[3*i+1], coord[3*i+2],
                  fmr->atom->fragment[i], // assuming next pivot state is state index 0, as in updatePivotState
                  veloc[3*i], veloc[3*i+1], veloc[3*i+2]
             );
    }
    fclose(fs);
  }
}


/*-----------------------------------------------------------------
  Subtracts translational and angular momentum from velocity 
-----------------------------------------------------------------*/
void Dynamics::removeTransRotFromVelocity()
{
  int natoms      = fmr->atom->natoms;
  double *coord   = fmr->atom->coord;
  double *veloc   = fmr->atom->veloc;
  double *mass    = fmr->atom->mass;
  int print_level = fmr->print_level;


  // Master does all the work, broadcasts later
  if (fmr->master_rank) {

    // ** Translate all coordinates to center of mass ** //
    printf("Translating center of mass to origin.\n");
    double COM[3];
    COM[0] = COM[1] = COM[2] = 0.0;
    for (int i=0; i<natoms; ++i) {
      COM[0] += mass[i] * coord[3*i]; 
      COM[1] += mass[i] * coord[3*i+1]; 
      COM[2] += mass[i] * coord[3*i+2]; 
    }
    COM[0] /= fmr->atom->totalMass;
    COM[1] /= fmr->atom->totalMass;
    COM[2] /= fmr->atom->totalMass;
    if (print_level > 0) {
      printf("Center of mass before:   %14.8f %14.8f %14.8f\n", COM[0], COM[1], COM[2]);
    }
    for (int i=0; i<natoms; ++i) {
      coord[3*i]   -= COM[0];
      coord[3*i+1] -= COM[1];
      coord[3*i+2] -= COM[2];
    }

    printf("Removing translation and rotation from velocities.\n");

    // ** Remove net translational momentum ** //
    double p[3];
    p[0] = p[1] = p[2] = 0.0;
    for (int i=0; i<natoms; ++i) {
      p[0] += mass[i] * veloc[3*i];    
      p[1] += mass[i] * veloc[3*i+1];    
      p[2] += mass[i] * veloc[3*i+2];    
    }
    if (print_level > 0) {
      printf("Net momentum before:     %14.8f %14.8f %14.8f\n", p[0], p[1], p[2]);
    }
    p[0] /= fmr->atom->totalMass;
    p[1] /= fmr->atom->totalMass;
    p[2] /= fmr->atom->totalMass;
    for (int i=0; i<natoms; ++i) {
      veloc[3*i]   -= p[0];
      veloc[3*i+1] -= p[1];
      veloc[3*i+2] -= p[2];
    }
    if (print_level > 0) {
      p[0] = p[1] = p[2] = 0.0;
      for (int i=0; i<natoms; ++i) {
        p[0] += mass[i] * veloc[3*i];    
        p[1] += mass[i] * veloc[3*i+1];    
        p[2] += mass[i] * veloc[3*i+2];    
      }
      printf("Net momentum after:      %14.8f %14.8f %14.8f\n", p[0], p[1], p[2]);
    }

    // ** Remove net angular momentum ** //
    // L = q x p
    double L[3];
    L[0] = L[1] = L[2] = 0.0;
    for (int i=0; i<natoms; ++i) {
      L[0] += mass[i] * (coord[3*i+1]*veloc[3*i+2] - coord[3*i+2]*veloc[3*i+1]);
      L[1] += mass[i] * (coord[3*i+2]*veloc[3*i  ] - coord[3*i  ]*veloc[3*i+2]);
      L[2] += mass[i] * (coord[3*i  ]*veloc[3*i+1] - coord[3*i+1]*veloc[3*i  ]);
    }
    if (print_level > 0) {
      printf("Angular momentum before: %14.8f %14.8f %14.8f\n", L[0], L[1], L[2]);
    }
    // Moment of inertia tensor I = m (q^2(diag) - q x q)
    double I[9];
    for (int i=0; i<9; ++i) I[i] = 0.0;
    for (int i=0; i<natoms; ++i) {
      const double x = coord[3*i];
      const double y = coord[3*i+1];
      const double z = coord[3*i+2];
      I[0*3 + 0] += mass[i] * (y*y + z*z);
      I[1*3 + 1] += mass[i] * (x*x + z*z);
      I[2*3 + 2] += mass[i] * (x*x + y*y);
      I[0*3 + 1] -= mass[i] * x*y;
      I[0*3 + 2] -= mass[i] * x*z;
      I[1*3 + 2] -= mass[i] * y*z;
    }
    I[1*3 + 0] = I[0*3 + 1];
    I[2*3 + 0] = I[0*3 + 2];
    I[2*3 + 1] = I[1*3 + 2];
    // Invert
    fmr->math->invertMatrix(I, 3);
    // Angular velocity, w = I^(-1).L
    double wx, wy, wz;
    wx = I[0*3+0]*L[0] + I[0*3+1]*L[1] + I[0*3+2]*L[2];
    wy = I[1*3+0]*L[0] + I[1*3+1]*L[1] + I[1*3+2]*L[2];
    wz = I[2*3+0]*L[0] + I[2*3+1]*L[1] + I[2*3+2]*L[2];
    if (print_level > 0) {
      printf("Angular velocity before: %14.8f %14.8f %14.8f\n", wx, wy, wz);
    }
    // corrected velocity v' = v - (w x q)
    for (int i=0; i<natoms; ++i) {
      const double x = coord[3*i];
      const double y = coord[3*i+1];
      const double z = coord[3*i+2];
      veloc[3*i  ] -= wy*z - wz*y;
      veloc[3*i+1] -= wz*x - wx*z;
      veloc[3*i+2] -= wx*y - wy*x;
    }
    if (print_level > 0) {
      // recompute the angular momentum, should be zero now
      L[0] = L[1] = L[2] = 0.0;
      for (int i=0; i<natoms; ++i) {
        L[0] += mass[i] * (coord[3*i+1]*veloc[3*i+2] - coord[3*i+2]*veloc[3*i+1]);
        L[1] += mass[i] * (coord[3*i+2]*veloc[3*i  ] - coord[3*i  ]*veloc[3*i+2]);
        L[2] += mass[i] * (coord[3*i  ]*veloc[3*i+1] - coord[3*i+1]*veloc[3*i  ]);
      }
      printf("Angular momentum after:  %14.8f %14.8f %14.8f\n", L[0], L[1], L[2]);
      wx = I[0*3+0]*L[0] + I[0*3+1]*L[1] + I[0*3+2]*L[2];
      wy = I[1*3+0]*L[0] + I[1*3+1]*L[1] + I[1*3+2]*L[2];
      wz = I[2*3+0]*L[0] + I[2*3+1]*L[1] + I[2*3+2]*L[2];
      printf("Angular velocity after:  %14.8f %14.8f %14.8f\n", wx, wy, wz);
    }

  }

  // ** Broadcast the adjusted coordinates and velocities ** //
  MPI_Bcast(coord, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
  MPI_Bcast(veloc, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);

}

/*-----------------------------------------------------------------
  Subtracts translational and angular momentum from gradient 
-----------------------------------------------------------------*/
void Dynamics::removeTransRotFromGradient()
{
  int natoms      = fmr->atom->natoms;
  double *coord   = fmr->atom->coord;
  double *mass    = fmr->atom->mass;
  double *force   = fmr->atom->force; // actually the gradient
  int print_level = fmr->print_level;

  // Master does all the work, broadcasts later
  if (fmr->master_rank) {
    
    // ** Remove COM translational force ** //
    double cx, cy, cz;
    cx = cy = cz = 0.0;
    for (int i=0; i<natoms; ++i) {
      cx += force[3*i]; 
      cy += force[3*i+1]; 
      cz += force[3*i+2]; 
    }
    cx /= fmr->atom->totalMass;
    cy /= fmr->atom->totalMass;
    cz /= fmr->atom->totalMass;
    for (int i=0; i<natoms; ++i) {
      force[3*i]   -= cx * mass[i];
      force[3*i+1] -= cy * mass[i];
      force[3*i+2] -= cz * mass[i];
    }
 
    // ** Compute net force now ** //

    // ** Remove net torque ** //
    // T = q x F  (mass removed)
    double tx, ty, tz;
    tx = ty = tz = 0.0;
    for (int i=0; i<natoms; ++i) {
      const double x  = coord[3*i];
      const double y  = coord[3*i+1];
      const double z  = coord[3*i+2];
      const double fx = force[3*i];
      const double fy = force[3*i+1];
      const double fz = force[3*i+2];
      tx += y*fz - z*fy;
      ty += z*fx - x*fz;
      tz += x*fy - y*fz;
    } 
    printf("Net torque before: %f %f %f\n", tx, ty, tz);
    // Need to remove individual atom contribution to torque
    // For a given atom, contribution is m_i.I^(-1), where I is inertia tensor
    // Moment of inertia tensor I = m (q^2(diag) - q x q)
    double I[9];
    for (int i=0; i<natoms; ++i) {
      const double x = coord[3*i];
      const double y = coord[3*i+1];
      const double z = coord[3*i+2];
      I[0*3 + 0] += mass[i] * (y*y + z*z);
      I[1*3 + 1] += mass[i] * (x*x + z*z);
      I[2*3 + 2] += mass[i] * (x*x + y*y);
      I[0*3 + 1] -= mass[i] * x*y;
      I[0*3 + 2] -= mass[i] * x*z;
      I[1*3 + 2] -= mass[i] * y*z;
    }
    I[1*3 + 0] = I[0*3 + 1];
    I[2*3 + 0] = I[0*3 + 2];
    I[2*3 + 1] = I[1*3 + 2];
    // Invert
    fmr->math->invertMatrix(I, 3);
    // Angular acceleration dw/dt = I^(-1).torque
    double Itx, Ity, Itz;
    Itx = I[0*3+0]*tx + I[0*3+1]*ty + I[0*3+2]*tz;
    Ity = I[1*3+0]*tx + I[1*3+1]*ty + I[1*3+2]*tz;
    Itz = I[2*3+0]*tx + I[2*3+1]*ty + I[2*3+2]*tz;
    // Compute correction to force: m_i * (dw/dt) x q
    for (int i=0; i<natoms; ++i) {
      const double x  = coord[3*i];
      const double y  = coord[3*i+1];
      const double z  = coord[3*i+2];
      force[3*i  ] -= mass[i] * (Ity*z - Itz*y);
      force[3*i+1] -= mass[i] * (Itz*x - Itx*z);
      force[3*i+2] -= mass[i] * (Itx*y - Ity*x);
    }

    // ** Printing if desired ** //
    if (print_level > 0) {

      // Center of mass: should remain at origin
      cx = cy = cz = 0.0;
      for (int i=0; i<natoms; ++i) {
        cx += mass[i] * coord[3*i]; 
        cy += mass[i] * coord[3*i+1]; 
        cz += mass[i] * coord[3*i+2]; 
      }
      cx /= fmr->atom->totalMass;
      cy /= fmr->atom->totalMass;
      cz /= fmr->atom->totalMass;

      // Net force
      double fx, fy, fz;
      fx = fy = fz = 0.0;
      for (int i=0; i<natoms; ++i) {
        fx += force[3*i];
        fy += force[3*i+1];
        fz += force[3*i+2];
      }

      // Net torque
      tx = ty = tz = 0.0;
      for (int i=0; i<natoms; ++i) {
        const double x  = coord[3*i];
        const double y  = coord[3*i+1];
        const double z  = coord[3*i+2];
        const double fx = force[3*i];
        const double fy = force[3*i+1];
        const double fz = force[3*i+2];
        tx += y*fz - z*fy;
        ty += z*fx - x*fz;
        tz += x*fy - y*fz;
      } 

      printf("Center of mass: %f %f %f\n", cx, cy, cz);
      printf("Net force:      %f %f %f\n", fx, fy, fz);
      printf("Net torque:     %f %f %f\n", tx, ty, tz);
      
    }
 
  }

  // ** Broadcast the adjusted force ** //
  MPI_Bcast(force, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);

}
