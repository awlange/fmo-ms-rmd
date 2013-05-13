/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "run.h"
#include "atom.h"
#include "matrix.h"
#include "dsyev.h"
#include "input.h"
#include "state.h"

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/

Matrix::Matrix(FMR *fmr) : Pointers(fmr)
{
  // Basic initializations
  GSEnergy   = 0.0;
  H          = NULL;
  Evecs      = NULL;
  Evals      = NULL;
  Repulsion  = NULL;
  GSCoeffs   = NULL;
  HX         = NULL;
  RepulsionX = NULL;
}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/

Matrix::~Matrix()
{
  int nstates = fmr->atom->nstates; 

  if (H != NULL) {
    for (int i=0; i<nstates; ++i)
      delete [] H[i];
    delete [] H; 
  }
  if (Evecs != NULL) {
    for (int i=0; i<nstates; ++i)
      delete [] Evecs[i];
    delete [] Evecs; 
  }
  if (Evals       != NULL) delete [] Evals; 
  if (Repulsion   != NULL) delete [] Repulsion; 
  if (GSCoeffs    != NULL) delete [] GSCoeffs; 
  if (HX != NULL) {
    for (int i=0; i<nstates; ++i) { 
      for (int j=0; j<nstates; ++j) {
        delete [] HX[i][j];
      }
      delete [] HX[i];
    } 
    delete [] HX; 
  }
  if (RepulsionX != NULL) {
    for (int i=0; i<nstates; ++i) {
      delete [] RepulsionX[i];
    }
    delete [] RepulsionX;
  }
}


/*-----------------------------------------------------------------
  Builds H using current FMO energies 
-----------------------------------------------------------------*/
void Matrix::buildH() {

  int nstates = fmr->atom->nstates; 
  int prev_nstates = fmr->atom->prev_nstates;

  // ** Allocate H matrix ** //
  // If H is already allocated from previous step, de-allocate
  if (H != NULL) {
    for (int i=0; i<prev_nstates; ++i) delete [] H[i];
    delete [] H; 
  }
  // Allocate H based on number of states for this step
  H = new double*[nstates];
  for (int i=0; i<nstates; ++i) {
    H[i] = new double[nstates];
  }

  // ** Initialize, put FMO energies into diagonal ** //
  for (int i=0; i<nstates; ++i) {
    H[i][i] = fmr->run->fmo_energies[i]; 
    for (int j=i+1; j<nstates; ++j) {
      H[i][j] = 0.0; 
      H[j][i] = 0.0; 
    }
  }


  // ** Compute the repulsion energy in each state ** // 
  // Allocate memory like for H
  if (Repulsion != NULL) delete [] Repulsion;
  Repulsion = new double[nstates];
  ComputeRepulsions();
  // Add to H diagonal
  for (int i=0; i<nstates; ++i) H[i][i] += Repulsion[i];

  // ** Compute the off-diagonal coupling ** //
  ComputeCouplings();

  if (fmr->print_level > 0) {
    if (fmr->master_rank) {
      printf("Hamiltonian matrix:\n");
      for (int i=0; i<nstates; ++i) {
        printf("%2d:", i);
        for (int j=0; j<nstates; ++j) {
          printf(" %16.8f", H[i][j]);
        }
        printf("\n");
      }
    }
  }
}


/*-----------------------------------------------------------------
  Diagonalizes H to get eigenvectors and eigenvalues 
  Note: H actually remains intact, unmodified, but that's fine
        b/c we only need the evecs and evals
-----------------------------------------------------------------*/
void Matrix::diagonalizeH() {

  int nstates = fmr->atom->nstates;
  int prev_nstates = fmr->atom->prev_nstates;

  // ** Allocate Evecs and Evals as needed ** //
  if (Evals != NULL) delete [] Evals;
  Evals = new double[nstates];   
  if (Evecs != NULL) {
    for (int i=0; i<prev_nstates; ++i) {
      delete [] Evecs[i];
    }
    delete [] Evecs;
  }
  Evecs = new double*[nstates];
  for (int i=0; i<nstates; ++i) {
    Evecs[i] = new double[nstates];
  }

  // ** Now call the handy diagonalization wrapper ** //
  if (fmr->master_rank) {
    printf("Diagonalizing H...\n");
    dsyev(H, nstates, Evals, Evecs);
  }
  MPI_Barrier(fmr->world);

  // ** Set which eigenvalue and eigenvector is the ground state, store ** //
  // Evals is sorted above such that minimum is element 0, Evecs corresponds
  GSEnergy = Evals[0];
  if (GSCoeffs != NULL) delete [] GSCoeffs; 
  GSCoeffs = new double[nstates];
  for (int i=0; i<nstates; ++i) GSCoeffs[i] = Evecs[i][0]; 

  // Master broadcasts results of diagonalization to all worker ranks
  MPI_Bcast(&GSEnergy,      1, MPI_DOUBLE, MASTER_RANK, fmr->world);
  MPI_Bcast(GSCoeffs, nstates, MPI_DOUBLE, MASTER_RANK, fmr->world);

  // ** Printing ** //
  if (fmr->master_rank) {
    // All state reporting
    if (fmr->print_level > 0) {
      printf("Eigenvalues:\n");
      for (int i=0; i<nstates; ++i) { 
        printf(" %d: %f hartree %f kcal/mol\n", i, Evals[i], Evals[i]*fmr->math->au2kcal);
      }
    }
    // Ground state reporting
    printf("FMO-MS-RMD energy: %f hartree = %f kcal/mol\n", GSEnergy, GSEnergy*fmr->math->au2kcal);
    if (fmr->print_level > 0) {
      printf("Ground state coefficients:\n");
      for (int i=0; i<nstates; ++i) {
        printf(" %f", GSCoeffs[i]);
      }
      printf("\n");
    }
  }

}


/*-----------------------------------------------------------------
  Builds HX using current FMO gradients 
-----------------------------------------------------------------*/
void Matrix::buildHX() {

  int nstates = fmr->atom->nstates; 
  int natoms  = fmr->atom->natoms;
  int prev_nstates = fmr->atom->prev_nstates;

  // ** Allocate HX matrix ** //
  // If HX is already allocated from previous step, de-allocate
  if (HX != NULL) {
    for (int i=0; i<prev_nstates; ++i) { 
      for (int j=0; j<prev_nstates; ++j) {
        delete [] HX[i][j];
      }
      delete [] HX[i];
    } 
    delete [] HX; 
  }
  // Allocate HX based on number of states for this step
  HX = new double**[nstates];
  for (int i=0; i<nstates; ++i) {
    HX[i] = new double*[nstates];
    for (int j=0; j<nstates; ++j) {
      HX[i][j] = new double[3*natoms];
    }
  }


  // ** Initialize, put FMO energies into diagonal ** //
  for (int I=0; I<nstates; ++I) {
    for (int i=0; i<natoms; ++i) {
      HX[I][I][3*i]   = fmr->run->fmo_gradients[I*3*natoms + 3*i]; 
      HX[I][I][3*i+1] = fmr->run->fmo_gradients[I*3*natoms + 3*i+1]; 
      HX[I][I][3*i+2] = fmr->run->fmo_gradients[I*3*natoms + 3*i+2]; 
    }
    for (int J=I+1; J<nstates; ++J) {
      for (int i=0; i<natoms; ++i) {
        HX[I][J][3*i] = HX[I][J][3*i+1] = HX[I][J][3*i+2] = 0.0; 
        HX[J][I][3*i] = HX[J][I][3*i+1] = HX[J][I][3*i+2] = 0.0; 
      }
    }
  }


  // ** Compute the repulsion energy in each state ** // 
  // Allocate memory like for H
  if (RepulsionX != NULL) {
    for (int I=0; I<prev_nstates; ++I) {
      delete [] RepulsionX[I];
    }
    delete [] RepulsionX;
  }
  RepulsionX = new double*[nstates];
  for (int I=0; I<nstates; ++I) {
    RepulsionX[I] = new double[3*natoms];
  }
  ComputeRepulsionsX();
  // Add to HX diagonal
  for (int I=0; I<nstates; ++I) {
    for (int i=0; i<natoms; ++i) { 
      HX[I][I][3*i]   += RepulsionX[I][3*i];
      HX[I][I][3*i+1] += RepulsionX[I][3*i+1];
      HX[I][I][3*i+2] += RepulsionX[I][3*i+2];
    }
  }

  // ** Compute the off-diagonal coupling ** //
  ComputeCouplingsX();

#ifdef FMR_DEBUG
  // ** Printing ** //
  if (fmr->print_level > 0) {
    if (fmr->master_rank) {
      printf("HX:\n");
      for (int I=0; I<nstates; ++I) {
        for (int J=0; J<nstates; ++J) {
          for (int i=0; i<natoms; ++i) {
            printf("H[%d][%d][%d]: %f %f %f\n", I, J, i, HX[I][J][3*i], HX[I][J][3*i+1], HX[I][J][3*i+2]);
          }
        }
      }
    }
  }
#endif

}


/*-----------------------------------------------------------------
  Computes the ground state gradient via Hellman-Feynman 
-----------------------------------------------------------------*/
void Matrix::ComputeHellmanFeynmanGradient() {

  int nstates = fmr->atom->nstates;
  int natoms  = fmr->atom->natoms;
  double *GSGradient = fmr->atom->force;  // Important: points to atom data

  // Only master rank computes
  if (fmr->master_rank) {

    for (int i=0; i<natoms; ++i) {
      GSGradient[3*i]   = 0.0;
      GSGradient[3*i+1] = 0.0;
      GSGradient[3*i+2] = 0.0;
    }
  
    // ** Using ground state coefficients, contract ** //
    for (int I=0; I<nstates; ++I) {
      double cI = GSCoeffs[I];
      for (int J=0; J<nstates; ++J) {
        double cJ = GSCoeffs[J];
        double cc = cI*cJ;
        for (int i=0; i<natoms; ++i) {
          GSGradient[3*i]   += cc * HX[I][J][3*i];
          GSGradient[3*i+1] += cc * HX[I][J][3*i+1];
          GSGradient[3*i+2] += cc * HX[I][J][3*i+2];
        }
      }
    }
  
    // ** Printing ** //
    if (fmr->print_level > 0) {
      printf("Ground state gradient:\n");
      for (int i=0; i<natoms; ++i) {
        printf("%3d %14.8f %14.8f %14.8f\n", i, GSGradient[3*i], GSGradient[3*i+1], GSGradient[3*i+2]);
      }
    }
  }

  // ** Broadcast the force to the other ranks ** //
  MPI_Bcast(GSGradient, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world); 
}
