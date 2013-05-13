/* AWGL */
#ifndef FMR_STATE_H
#define FMR_STATE_H

#include "pointers.h"

namespace FMR_NS {

class State : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   State(FMR *);
   ~State();

   // ** Variables ** //
   int next_pivot_state; // The state index of the next step's pivot state
   int max_hops;         // Maximum number of hops in search
   double cut_OH;        // Distance cutoff in state search between O and H atoms
   int flag_read_MOs;    // Flag to indicate if it is safe to read MO coefficients from file from previous step
   int flag_state_number_change; // Flag to indicate a change in the number of states b/w steps

   // ** Functions ** //
   void state_search();           // The breadth-first search for fragmentation states
   void write_qchem_inputs(int);  // Writes the inputs for Q-Chem
   void updatePivotState();       // Update pivot state information *after* matrix diagonalization

};

}

#endif
