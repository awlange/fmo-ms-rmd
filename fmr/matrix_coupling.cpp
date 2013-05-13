/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "matrix.h"
#include "math.h"

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Compute the coupling energies, store in H matrix 
-----------------------------------------------------------------*/
void Matrix::ComputeCouplings() {
  // For now, only the master rank computes off-diagonal
  // Note: worker ranks do not have coupling stored in matrix for now
  if (fmr->master_rank) {
    printf("Computing coupling...\n");
    for (int I=0; I<fmr->atom->nstates; ++I) {
      for (int J=I+1; J<fmr->atom->nstates; ++J) {
        H[I][J] = ComputeCoupling(I,J);
        H[J][I] = H[I][J];
      }
    }
  }
}

/*-----------------------------------------------------------------
  Compute the coupling energy b/w states I and J 
-----------------------------------------------------------------*/
double Matrix::ComputeCoupling(int I, int J) {


  // Make things shorter hand here
  int natoms = fmr->atom->natoms;
  double *xx = fmr->atom->coord;
  char *sym  = fmr->atom->symbol;
  int *react = fmr->atom->reactive;
  double cut_OH  = fmr->state->cut_OH;
  double cut_OH2 = fmr->state->cut_OH * fmr->state->cut_OH;
  double inv_cut_OH6 = 1.0 / (cut_OH2 * cut_OH2 * cut_OH2);

  double ecouple = 0.0;

  // If not shared atom b/w I and J, then coupling is zero 
  int shared = 0;
  int iproton = -1; // index of shuttling proton b/w states I and J
  for (int i=0; i<natoms; ++i) {
    if (react[I*natoms + i] && react[J*natoms + i] && sym[i] == 'H') {
      // Reactive fragment in I shares at least one H atom with reactive fragment in J
      shared = 1;
      iproton = i;
      break;
    } 
  } 
  if (shared) {
    // Must check that I and J reactive fragments are not on the same oxygen,
    // as it is not possible for self-coupling
    // This can occur with bifurcated waters
    for (int i=0; i<natoms; ++i) {
      if (react[I*natoms + i] && react[J*natoms + i] && sym[i] == 'O') {
        // Same oxygen! No coupling allowed.
        shared = 0;
        break;
      } 
    }
  }

  // ** If passed above, now actually compute the coupling ** //
  if (shared) {
    // Compute interstate coupling repulsion
    double erepItoJ = 0.0;
    for (int i=0; i<natoms; ++i) {
      if (sym[i] == 'H' && react[I*natoms + i] && i == iproton) { // reactive H in state I
        for (int j=0; j<natoms; ++j) {
          if (sym[j] == 'O' && react[J*natoms + j]) { // reactive O in state J
            double dx = xx[3*i]   - xx[3*j];
            double dy = xx[3*i+1] - xx[3*j+1];
            double dz = xx[3*i+2] - xx[3*j+2];
            double dd = dx*dx + dy*dy + dz*dz;
            if (dd < cut_OH2) {
              double tmp = Sr(dd, 0.0, cut_OH2); 
              //double tmp = cut_OH / sqrt(dd) - 1.0;
              erepItoJ += tmp*tmp; 
            } 
          }
        }
      }
    }
    double erepJtoI = 0.0;
    for (int i=0; i<natoms; ++i) {
      if (sym[i] == 'H' && react[J*natoms + i] && i == iproton) { // reactive H in state J
        for (int j=0; j<natoms; ++j) {
          if (sym[j] == 'O' && react[I*natoms + j]) { // reactive O in state I
            double dx = xx[3*i]   - xx[3*j];
            double dy = xx[3*i+1] - xx[3*j+1];
            double dz = xx[3*i+2] - xx[3*j+2];
            double dd = dx*dx + dy*dy + dz*dz;
            if (dd < cut_OH2) {
              double tmp = Sr(dd, 0.0, cut_OH2); 
              //double tmp = cut_OH / sqrt(dd) - 1.0;
              erepJtoI += tmp*tmp; 
            } 
          }
        }
      }
    }
    // Geometric mean 
    ecouple = fmr->math->AA * sqrt(erepItoJ * erepJtoI);
  }

  return ecouple;
}


/*-----------------------------------------------------------------
  Compute the coupling gradients, store in HX matrix 
-----------------------------------------------------------------*/
void Matrix::ComputeCouplingsX() {
  // TODO
  // For now, only the master rank computes off-diagonal
  // Note: worker ranks do not have coupling stored in matrix for now
  if (fmr->master_rank) {
    printf("Computing coupling gradient...\n");
    for (int I=0; I<fmr->atom->nstates; ++I) {
      for (int J=I+1; J<fmr->atom->nstates; ++J) {
        ComputeCouplingX(I,J);
      }
    }
  }
}

/*-----------------------------------------------------------------
  Compute the coupling energy b/w states I and J 
-----------------------------------------------------------------*/
void Matrix::ComputeCouplingX(int I, int J) {

  // Make things shorter hand here
  int natoms = fmr->atom->natoms;
  double *xx = fmr->atom->coord;
  char *sym  = fmr->atom->symbol;
  int *react = fmr->atom->reactive;
  double cut_OH  = fmr->state->cut_OH;
  double cut_OH2 = fmr->state->cut_OH * fmr->state->cut_OH;
  double inv_cut_OH6 = 1.0 / (cut_OH2 * cut_OH2 * cut_OH2);

  // If not shared atom b/w I and J, then coupling is zero 
  int shared = 0;
  int iproton = -1; // index of shuttling proton b/w states I and J
  for (int i=0; i<natoms; ++i) {
    if (react[I*natoms + i] && react[J*natoms + i] && sym[i] == 'H') {
      // Reactive fragment in I shares at least one H atom with reactive fragment in J
      shared = 1;
      iproton = i;
      break;
    } 
  } 
  if (shared) {
    // Must check that I and J reactive fragments are not on the same oxygen,
    // as it is not possible for self-coupling
    // This can occur with bifurcated waters
    for (int i=0; i<natoms; ++i) {
      if (react[I*natoms + i] && react[J*natoms + i] && sym[i] == 'O') {
        // Same oxygen! No coupling allowed.
        shared = 0;
        break;
      } 
    }
  }

  if (shared) {
    // Compute interstate coupling repulsion, shows up as a chain rule factor
    double erepItoJ = 0.0;
    for (int i=0; i<natoms; ++i) {
      if (sym[i] == 'H' && react[I*natoms + i] && i == iproton) { // reactive H in state I
        for (int j=0; j<natoms; ++j) {
          if (sym[j] == 'O' && react[J*natoms + j]) { // reactive O in state J
            double dx = xx[3*i]   - xx[3*j];
            double dy = xx[3*i+1] - xx[3*j+1];
            double dz = xx[3*i+2] - xx[3*j+2];
            double dd = dx*dx + dy*dy + dz*dz;
            if (dd < cut_OH2) {
              double tmp = Sr(dd, 0.0, cut_OH2); 
              //double tmp = cut_OH / sqrt(dd) - 1.0;
              erepItoJ += tmp*tmp; 
            } 
          }
        }
      }
    }
    double erepJtoI = 0.0;
    for (int i=0; i<natoms; ++i) {
      if (sym[i] == 'H' && react[J*natoms + i] && i == iproton) { // reactive H in state J
        for (int j=0; j<natoms; ++j) {
          if (sym[j] == 'O' && react[I*natoms + j]) { // reactive O in state I
            double dx = xx[3*i]   - xx[3*j];
            double dy = xx[3*i+1] - xx[3*j+1];
            double dz = xx[3*i+2] - xx[3*j+2];
            double dd = dx*dx + dy*dy + dz*dz;
            if (dd < cut_OH2) {
              double tmp = Sr(dd, 0.0, cut_OH2); 
              //double tmp = cut_OH / sqrt(dd) - 1.0;
              erepJtoI += tmp*tmp; 
            } 
          }
        }
      }
    }
    double prefactor = fmr->math->AA / (2.0 * sqrt(erepItoJ * erepJtoI)); 

    // Now get the gradients
    for (int i=0; i<natoms; ++i) {
      if (sym[i] == 'H' && react[I*natoms + i] && i == iproton) { // reactive H in state I
        for (int j=0; j<natoms; ++j) {
          if (sym[j] == 'O' && react[J*natoms + j]) { // reactive O in state J
            double dx = xx[3*i]   - xx[3*j];
            double dy = xx[3*i+1] - xx[3*j+1];
            double dz = xx[3*i+2] - xx[3*j+2];
            double dd = dx*dx + dy*dy + dz*dz;
            if (dd < cut_OH2) {
              double d = sqrt(dd);
              double tmp = 2.0 * Sr(dd, 0.0, cut_OH2) * dSrdr(dd, 0.0, cut_OH2);
              //double tmp = 2.0 * (cut_OH/d - 1.0) * (-cut_OH/dd);
              tmp = prefactor * erepJtoI * tmp; 
              // Convert from hartree/Angs to hartree/bohr
              tmp *= fmr->math->b2a;
              double fx = tmp * dx / d;
              double fy = tmp * dy / d;
              double fz = tmp * dz / d;
              HX[I][J][3*i]   += fx;
              HX[I][J][3*i+1] += fy;
              HX[I][J][3*i+2] += fz;
              HX[I][J][3*j]   -= fx;
              HX[I][J][3*j+1] -= fy;
              HX[I][J][3*j+2] -= fz;
            } 
          }
        }
      }
    }
    for (int i=0; i<natoms; ++i) {
      if (sym[i] == 'H' && react[J*natoms + i]  && i == iproton) { // reactive H in state J
        for (int j=0; j<natoms; ++j) {
          if (sym[j] == 'O' && react[I*natoms + j]) { // reactive O in state I
            double dx = xx[3*i]   - xx[3*j];
            double dy = xx[3*i+1] - xx[3*j+1];
            double dz = xx[3*i+2] - xx[3*j+2];
            double dd = dx*dx + dy*dy + dz*dz;
            if (dd < cut_OH2) {
              double d = sqrt(dd);
              double tmp = 2.0 * Sr(dd, 0.0, cut_OH2) * dSrdr(dd, 0.0, cut_OH2);
              //double tmp = 2.0 * (cut_OH/d - 1.0) * (-cut_OH/dd);
              tmp = prefactor * erepItoJ * tmp; 
              // Convert from hartree/Angs to hartree/bohr
              tmp *= fmr->math->b2a;
              double fx = tmp * dx / d;
              double fy = tmp * dy / d;
              double fz = tmp * dz / d;
              HX[I][J][3*i]   += fx;
              HX[I][J][3*i+1] += fy;
              HX[I][J][3*i+2] += fz;
              HX[I][J][3*j]   -= fx;
              HX[I][J][3*j+1] -= fy;
              HX[I][J][3*j+2] -= fz;
            } 
          }
        }
      }
    }

    // ** Symmetrize coupling gradient for convenience ** //
    for (int i=0; i<natoms; ++i) {
      HX[J][I][3*i]   = HX[I][J][3*i];
      HX[J][I][3*i+1] = HX[I][J][3*i+1];
      HX[J][I][3*i+2] = HX[I][J][3*i+2];
    }
  }

}
