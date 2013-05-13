/* AWGL */
#ifndef FMR_RUN_H
#define FMR_RUN_H

#include "pointers.h"

// Run type definitions
#define RUN_ENERGY   1
#define RUN_FORCE    2
#define RUN_MOLDYN   3
#define RUN_TYPE_MIN 1
#define RUN_TYPE_MAX 3

namespace FMR_NS {


class Run : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Run(FMR *);
   ~Run();

   // ** Variables ** //
   int run_type;              // Run type index
   int FMO_only;              // Run with only FMO, no Multi-state FMO
   int EnvApprox;             // Environment approximation on/off

   char qchem_scratch[256];   // Full path of Q-Chem scratch directory 
   char qchem_exec[256];      // Full path of Q-Chem exectuable

   int n_monomers;            // (# fragments) * (# states)
   int n_dimers;              // ((# fragments) * (# fragments - 1) / 2 ) * (# states)
   int n_dimers_sq;           // ((# fragments) * (# fragments) ) * (# states)
   double *fmo_energies;      // FMO energies for each state
   double *monomer_energies;  // FMO monomer energies
   double *dimer_energies;    // FMO dimer energies
   double *fmo_gradients;     // FMO gradients for each state
   double *monomer_gradients; // FMO monomer gradients
   double *dimer_gradients;   // FMO dimer gradients

   // ** Functions ** //
   void run_calculation();     // General run calculation
   void calculate_energy();    // Computes energy only
   void calculate_force();     // Computes energy + force
   void calculate_moldyn();    // Computes molecular dynamics
   void do_fmo_calculations(int); // Performs the all FMO calculations in parallel
   void do_fmo_calculations_env();  // Performs the FMO calculations in parallel, using env approximation
   void perturb_coords();      // Just to perturb the coordinates slightly randomly

};

}

#endif
