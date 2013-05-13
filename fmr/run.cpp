/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "run.h"
#include "atom.h"
#include "state.h"
#include "input.h"
#include "matrix.h"
#include "dynamics.h"

#define MAX_LENGTH 1024

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/

Run::Run(FMR *fmr) : Pointers(fmr)
{
   // Set defaults
   run_type   = RUN_ENERGY;
   FMO_only   = 0; 
   EnvApprox  = 0;
   sprintf(qchem_scratch, "%s", "null");
   sprintf(qchem_exec, "%s", "null");
 
   n_monomers = 0;
   n_dimers   = 0;
   fmo_energies      = NULL;
   monomer_energies  = NULL;
   dimer_energies    = NULL;
   fmo_gradients     = NULL;
   monomer_gradients = NULL;
   dimer_gradients   = NULL;
}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/

Run::~Run()
{
  if (fmo_energies      != NULL) delete [] fmo_energies;
  if (monomer_energies  != NULL) delete [] monomer_energies;
  if (dimer_energies    != NULL) delete [] dimer_energies;
  if (fmo_gradients     != NULL) delete [] fmo_gradients;
  if (monomer_gradients != NULL) delete [] monomer_gradients;
  if (dimer_gradients   != NULL) delete [] dimer_gradients;
}

/*-----------------------------------------------------------------
  Run calculation: Upper level run command 
-----------------------------------------------------------------*/
void Run::run_calculation()
{
  // ** Start clock ** //
  MPI_Barrier(fmr->world);
  double clock_start = MPI_Wtime();

  // ** Do the calculation ** // 
  if (run_type == RUN_ENERGY) {
    calculate_energy();
  }
  else if (run_type == RUN_FORCE) {
    calculate_force();
  }
  else if (run_type == RUN_MOLDYN) {
    calculate_moldyn();
  }


  // ** Stop clock ** //
  MPI_Barrier(fmr->world);
  double clock_end = MPI_Wtime();
  double run_time = clock_end - clock_start;
  if (fmr->master_rank) {
    printf("\nTotal run time: %20.6f seconds on %d ranks\n", run_time, fmr->world_size);
  }
}


/*-----------------------------------------------------------------
  Calculate energy 
-----------------------------------------------------------------*/
void Run::calculate_energy()
{
  // Step 1. Do state search to determine fragmentation states
  fmr->state->state_search();  
  // Step 2. Write Q-Chem inputs for all FMO calculations
  fmr->state->write_qchem_inputs(RUN_ENERGY);
  // Step 3. Divide up FMO calculation and run in parallel
  do_fmo_calculations(0);
  if (!fmr->run->FMO_only) {
    // Step 4. Construct model Hamiltonian
    fmr->matrix->buildH();
    // Step 5. Diagonalize Hamiltonian
    fmr->matrix->diagonalizeH();
  }
}

/*-----------------------------------------------------------------
  Calculate force (well, really the gradient)
-----------------------------------------------------------------*/
void Run::calculate_force()
{
  // Step 1. Do state search to determine fragmentation states
  fmr->state->state_search();  
  // Step 2. Write Q-Chem inputs for all FMO calculations
  fmr->state->write_qchem_inputs(RUN_FORCE);
  // Step 3. Divide up FMO calculation and run in parallel
  do_fmo_calculations(1);
  // Step 4. Construct model Hamiltonian
  if (!fmr->run->FMO_only) {
    // Step 4. Construct model Hamiltonian
    fmr->matrix->buildH();
    // Step 5. Diagonalize Hamiltonian
    fmr->matrix->diagonalizeH();
    // Step 6. Build HX for Hellman-Feynman 
    fmr->matrix->buildHX();
    // Step 7. Compute ground state force via Hellman-Feynman
    fmr->matrix->ComputeHellmanFeynmanGradient();
    // Step 8. Update pivot state information, etc. for next step
    fmr->state->updatePivotState();
  } else {
    // Well, I haven't put this in here yet, but you can do it by just turning off the repulsion and coupling
    if (fmr->master_rank) printf("Sorry. FMO-only gradient not yet in place.\n");
    fmr->error(FLERR, "Sorry. FMO-only gradient not yet in place.\n");
  }
}

/*-----------------------------------------------------------------
  Calculate molecular dynamics 
-----------------------------------------------------------------*/
void Run::calculate_moldyn()
{
  // Step 1. Calculate force
  calculate_force(); 
  // Step 2. Initialize for MD
  fmr->dynamics->init();
  // Step 3. Run MD for number of steps
  fmr->dynamics->runMolecularDynamics();
}


/*-----------------------------------------------------------------
  Perturb coordinates 
-----------------------------------------------------------------*/
void Run::perturb_coords()
{
  // Seed
  fmr->math->rng_seed(999); // same for now
  double max_dist = 0.0001; // Angs

  int natoms = fmr->atom->natoms;
  for (int i=0; i<natoms; ++i) {
     // add plus/minus max_dist
     atom->coord[3*i]   += max_dist * (2.0*fmr->math->rng() - 1.0); 
     atom->coord[3*i+1] += max_dist * (2.0*fmr->math->rng() - 1.0); 
     atom->coord[3*i+2] += max_dist * (2.0*fmr->math->rng() - 1.0);
  }
}
