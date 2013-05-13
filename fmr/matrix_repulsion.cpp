/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "matrix.h"
#include "math.h"

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Compute the repulsion energy for all states, store in memory 
-----------------------------------------------------------------*/
void Matrix::ComputeRepulsions() {
  // For now, just make the master rank do this
  if (fmr->master_rank) {
    printf("Computing repulsion...\n");
    for (int I=0; I<fmr->atom->nstates; ++I) {
      Repulsion[I] =  ComputeRepulsionForState(I);
      Repulsion[I] += ComputeInternalPenaltyForState(I);
    }
    if (fmr->print_level > 0) {
      printf("Repulsion energies:\n");
      for (int I=0; I<fmr->atom->nstates; ++I) {
        printf("%2d: %16.8f\n", I, Repulsion[I]);
      }
    }
  }
  //MPI_Bcast(fmr->matrix->Repulsion, fmr->atom->nstates, MPI_DOUBLE, MASTER_RANK, fmr->world);
}

/*-----------------------------------------------------------------
  Compute the repulsion energy for a specified state I 
-----------------------------------------------------------------*/
double Matrix::ComputeRepulsionForState(int I) {

  // Make things shorter hand here
  int natoms = fmr->atom->natoms;
  double *xx = fmr->atom->coord;
  char *sym  = fmr->atom->symbol;
  int *react = fmr->atom->reactive;
  double cut_OH  = fmr->state->cut_OH;
  double cut_OH2 = fmr->state->cut_OH * fmr->state->cut_OH;

  double erep = 0.0;

  // Let each hydronium H atom interact withh all non-hydronium O atoms in state I
  for (int i=0; i<natoms; ++i) {
    if (sym[i] == 'H' && react[I*natoms + i]) {
      for (int j=0; j<natoms; ++j) {
        if (sym[j] == 'O' && react[I*natoms + j] == 0) {
          double dx = xx[3*i]   - xx[3*j];
          double dy = xx[3*i+1] - xx[3*j+1];
          double dz = xx[3*i+2] - xx[3*j+2];
          double dd = dx*dx + dy*dy + dz*dz;
          if (dd < cut_OH2) {
            double tmp = Sr(dd, 0.0, cut_OH2);
            //double tmp = RR / sqrt(dd) - 1.0;
            erep += tmp*tmp; 
          }
        }
      }
    }
  }

  return erep * fmr->math->AA;
}

/*-----------------------------------------------------------------
  Compute the repulsion gradients for all states, store in memory 
-----------------------------------------------------------------*/
void Matrix::ComputeRepulsionsX() {
  // For now, just make the master rank do this
  if (fmr->master_rank) {
    printf("Computing repulsion gradients...\n");
    for (int I=0; I < fmr->atom->nstates; ++I) {
      // Clear out for this state
      for (int i=0; i < fmr->atom->natoms; ++i) {
        RepulsionX[I][3*i]   = 0.0;
        RepulsionX[I][3*i+1] = 0.0;
        RepulsionX[I][3*i+2] = 0.0;
      }
      ComputeRepulsionXForState(I);
      ComputeInternalPenaltyXForState(I);
    }
  }
  //MPI_Bcast(fmr->matrix->Repulsion, fmr->atom->nstates, MPI_DOUBLE, MASTER_RANK, fmr->world);
}

/*-----------------------------------------------------------------
  Compute the repulsion energy for a specified state I 
-----------------------------------------------------------------*/
void Matrix::ComputeRepulsionXForState(int I) {

  // Make things shorter hand here
  int natoms = fmr->atom->natoms;
  double *xx = fmr->atom->coord;
  char *sym  = fmr->atom->symbol;
  int *react = fmr->atom->reactive;
  double cut_OH  = fmr->state->cut_OH;
  double cut_OH2 = fmr->state->cut_OH * fmr->state->cut_OH;

  // Let each hydronium H atom interact withh all non-hydronium O atoms in state I
  for (int i=0; i<natoms; ++i) {
    if (sym[i] == 'H' && react[I*natoms + i]) {
      for (int j=0; j<natoms; ++j) {
        if (sym[j] == 'O' && react[I*natoms + j] == 0) {
          double dx = xx[3*i]   - xx[3*j];
          double dy = xx[3*i+1] - xx[3*j+1];
          double dz = xx[3*i+2] - xx[3*j+2];
          double dd = dx*dx + dy*dy + dz*dz;
          if (dd < cut_OH2) {
            double d = sqrt(dd);
            double tmp = 2.0 * Sr(dd, 0.0, cut_OH2) * dSrdr(dd, 0.0, cut_OH2); 
            //double tmp = 2.0 * (RR/d - 1.0) * (-RR/dd); 
            tmp = fmr->math->AA * tmp; 
            // Convert from hartree/Angs to hartree/bohr
            tmp *= fmr->math->b2a;
            double fx = tmp * dx / d;
            double fy = tmp * dy / d;
            double fz = tmp * dz / d;
            RepulsionX[I][3*i]   += fx;
            RepulsionX[I][3*i+1] += fy;
            RepulsionX[I][3*i+2] += fz;
            RepulsionX[I][3*j]   -= fx;
            RepulsionX[I][3*j+1] -= fy;
            RepulsionX[I][3*j+2] -= fz;
          }
        }
      }
    }
  }

}


/*-----------------------------------------------------------------
  Compute the internal penalty energy for a specified state I 
-----------------------------------------------------------------*/
double Matrix::ComputeInternalPenaltyForState(int I) {

  // Make things shorter hand here
  int natoms = fmr->atom->natoms;
  double *xx = fmr->atom->coord;
  char *sym  = fmr->atom->symbol;
  int *react = fmr->atom->reactive;

  double BB = fmr->math->BB; 
  double r0 = fmr->math->rb;
  double r02 = r0*r0;
  double cut_OH2 = fmr->state->cut_OH * fmr->state->cut_OH;
  double denom = cut_OH2 - r02;
  denom *= denom * denom;

  double eint = 0.0;

  for (int i=0; i<natoms; ++i) {
    if (sym[i] == 'H') {
      for (int j=0; j<natoms; ++j) {
        if (sym[j] == 'O') {
          if (fmr->atom->fragment[I*natoms + i] == fmr->atom->fragment[I*natoms + j]) { // same fragment?
            double dx = xx[3*i]   - xx[3*j];
            double dy = xx[3*i+1] - xx[3*j+1];
            double dz = xx[3*i+2] - xx[3*j+2];
            double dd = dx*dx + dy*dy + dz*dz;
            if (dd > r02) {
              double tmp = 1.0 - Sr(dd, r02, cut_OH2); 
              eint += tmp*tmp; 
            }
          }
        }
      }
    }
  }

  printf("Penalty for state %d: %f\n", I, eint*BB);

  return eint * BB; 
}


/*-----------------------------------------------------------------
  Compute the internal penalty gradient for a specified state I 
-----------------------------------------------------------------*/
void Matrix::ComputeInternalPenaltyXForState(int I) {

  // Make things shorter hand here
  int natoms = fmr->atom->natoms;
  double *xx = fmr->atom->coord;
  char *sym  = fmr->atom->symbol;
  int *react = fmr->atom->reactive;

  double BB = fmr->math->BB; 
  double r0 = fmr->math->rb;
  double r02 = r0*r0;
  double cut_OH2 = fmr->state->cut_OH * fmr->state->cut_OH;
  double denom = cut_OH2 - r02;
  denom *= denom * denom;

  for (int i=0; i<natoms; ++i) {
    if (sym[i] == 'H') {
      for (int j=0; j<natoms; ++j) {
        if (sym[j] == 'O') {
          if (fmr->atom->fragment[I*natoms + i] == fmr->atom->fragment[I*natoms + j]) { // same fragment?
            double dx = xx[3*i]   - xx[3*j];
            double dy = xx[3*i+1] - xx[3*j+1];
            double dz = xx[3*i+2] - xx[3*j+2];
            double dd = dx*dx + dy*dy + dz*dz;
            if (dd > r02) {
              double tmp = 2.0 * (1.0 - Sr(dd, r02, cut_OH2)) * (-dSrdr(dd, r02, cut_OH2));
              tmp = BB * tmp; 
              // Convert from hartree/Angs to hartree/bohr
              tmp *= fmr->math->b2a;
              double d = sqrt(dd);
              double fx = tmp * dx / d;
              double fy = tmp * dy / d;
              double fz = tmp * dz / d;
              RepulsionX[I][3*i]   += fx;
              RepulsionX[I][3*i+1] += fy;
              RepulsionX[I][3*i+2] += fz;
              RepulsionX[I][3*j]   -= fx;
              RepulsionX[I][3*j+1] -= fy;
              RepulsionX[I][3*j+2] -= fz;
            }
          }
        }
      }
    }
  }
}
