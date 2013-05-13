/* AWGL */
#include "fmr.h"
#include "math.h"
#include "fmr_math.h"
#include "stdlib.h"
#include "time.h"

// ** LAPACK functions used for matrix inversion ** //
extern "C" {
  // LU decomoposition of a general matrix
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  // generate inverse of a matrix given its LU decomposition
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/
Math::Math(FMR *fmr) : Pointers(fmr)
{
  // Math constants 
  pi = 3.1415926535897932384; 
  twopi = 2.0 * pi;

  // Converson factors
  b2a     = 0.529177249;
  a2b     = 1.0/b2a;
  au2kcal = 627.50947;
  kcal2au = 1.0/au2kcal;
  kB_eV   = 0.00008617332478;
  kB_au   = kB_eV / 27.211396132;
  amu2au  = 1822.86045493;
  fem2au  = 1.0 / 0.02418884326505;
  seed    = 999; // default
  hasBeenSeeded = false;

  // FMO-MS-RMD parameters
  // cc-pVDZ/MP2
  //AA = 0.0083290406; // hartree    // Geometric mean, force continuity 
  AA = 0.0064054303;
  BB = 0.8; // hartree / Angs**2 
  CC = 0.0; 
  rb = 1.4; // Angs

}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/
Math::~Math()
{

}

/* ----------------------------------------------------------------------------
  Simple random number generator seed 
---------------------------------------------------------------------------- */

void Math::rng_seed (int n) 
{
  if (!hasBeenSeeded) {
    int n_used = n;
    if (n < 0) n_used = time(NULL);
    srandom(n_used);
    if (fmr->master_rank) {
      printf("Random number generator seeded with: %d\n", n_used);
    }
    hasBeenSeeded = true;
  }
}

/* ----------------------------------------------------------------------------
  Simple random number generator 
  Produces a number in range [0:1]
---------------------------------------------------------------------------- */

double Math::rng() 
{
  return random() / (double)RAND_MAX;
}

/* ----------------------------------------------------------------------------
  Produces a number in range [min:max]
---------------------------------------------------------------------------- */
double Math::rngInRange(double min, double max)
{
  return min + rng() * (max-min);
}

/* ----------------------------------------------------------------------------
  Inverts matrix A which is NxN in size 
  After call, A is inverted
  This was adapted from a posting on stackoverflow.com
---------------------------------------------------------------------------- */
void Math::invertMatrix(double *A, int N)
{
  int *IPIV = new int[N+1];
  int LWORK = N*N;
  double *WORK = new double[LWORK];
  int INFO;

  dgetrf_(&N,&N,A,&N,IPIV,&INFO);
  dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

  delete IPIV;
  delete WORK;
}
