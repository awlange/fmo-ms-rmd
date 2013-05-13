/* AWGL */
#include "fmr.h"

using namespace FMR_NS;

/* ----------------------------------------------------------------------------
  The main program
---------------------------------------------------------------------------- */

int main(int argc, char **argv) 
{

  MPI_Init(&argc, &argv);
  
  FMR *fmr = new FMR(argc, argv, MPI_COMM_WORLD);
  fmr->execute();
  delete fmr;

  MPI_Finalize();
  return 0;
}
