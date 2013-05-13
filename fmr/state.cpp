/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "run.h"
#include "matrix.h"
#include <vector>

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/
State::State(FMR *fmr) : Pointers(fmr)
{
  // Basic initializations
  next_pivot_state = 0;
  //cut_OH   = 2.5; // Angstroms
  cut_OH   = 2.4; // Angstroms
  max_hops = 2;
  flag_read_MOs = 0; // default to not okay to read
  flag_state_number_change = 0;
}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/
State::~State()
{

}

/*-----------------------------------------------------------------
  State search algorithm 
-----------------------------------------------------------------*/
void State::state_search()
{
  // Given the atoms and an initial (pivot) fragmentation,
  // locate other possible states 

  int natoms = fmr->atom->natoms;
  int fragmentIndexChange = 0;

  // Previous step's number of states
  fmr->atom->prev_nstates = fmr->atom->nstates; 

  if (fmr->master_rank) {
    printf("Beginning state search on rank %d\n", fmr->my_rank);

    Atom *atom     = fmr->atom;
    int nfragments = atom->nfragments;
    int ireactive  = atom->ireactive; // which fragment is reactive donor
    int nstates    = 1; // # states found
    int state      = 0; // Current state index from which we are searching
    int nextstate  = 1; // Next state index in search
    int hop        = 0; // Hop level
    int hop_length = 1; // How many hop elements we've counted
    double cut_OH2 = cut_OH * cut_OH;

    // Initialize availability and reactive fragment
    // Here, we assume that state zero is the pivot state
    for (int i=0; i<natoms; ++i) {
      atom->available[i] = 1; // Default all to be available
      atom->hop[i] = 0;       // Default to zero, pivot state 
      int ifrag = atom->fragment[i];
      if (ifrag == ireactive) {
	atom->reactive[i] = 1;
	if (atom->symbol[i] == 'O') {
	  // Hydronium oxygen is unavailable to accept
	  atom->available[i] = 0;
	}
      } else {
	atom->reactive[i] = 0;
      }
    }

    // Helper vector for keeping track of which oxygens are unavailable
    std::vector<int> toBecomeUnavailableNextHop;
    // Helper array for checking that fragment indexing has changed or not
    int *old_fragment = new int [natoms*MAX_STATES];
    for (int i=0; i<natoms*MAX_STATES; ++i) {
      old_fragment[i] = atom->fragment[i];
    }

    // *** Breadth first search *** //
    // Loop until one of following is met:
    // 1. Exhausted all possible donations
    // 2. Maximum hops reached
    // 3. Maximum states reached
    while (state < nextstate && state < MAX_STATES && hop < max_hops) {


      // Loop through atoms of state and attempt to donate
      // hydronium H's to all possible O's within cutoff that are available
      for (int i=0; i<natoms; ++i) {
	if (atom->symbol[i] == 'H' && atom->available[i] && atom->reactive[state*natoms + i]) {
	  for (int j=0; j<natoms; ++j) {
	    if (atom->symbol[j] == 'O' && atom->available[j] && !atom->reactive[state*natoms + j]) {
	      double dx = atom->coord[3*i]   - atom->coord[3*j];
	      double dy = atom->coord[3*i+1] - atom->coord[3*j+1];
	      double dz = atom->coord[3*i+2] - atom->coord[3*j+2];
	      double dd = dx*dx + dy*dy + dz*dz;
	      if (dd < cut_OH2) {
		// Copy current state to next state, setting acceptor fragment
		// as the hydronium fragment in the next state
		for (int k=0; k<natoms; ++k) {
		  atom->fragment[nstates*natoms + k] = atom->fragment[state*natoms + k];
		  atom->hop[nstates*natoms + k]      = atom->hop[state*natoms + k] + 1; // add one hop
                  hop_length++;
		  if (atom->fragment[state*natoms + k] == atom->fragment[state*natoms + j]) {
		    atom->reactive[nstates*natoms + k] = 1;
		  } else {
		    atom->reactive[nstates*natoms + k] = 0;
		  }
		}
		// Swap fragment index of donated H atom and set to hydronium
		atom->fragment[nextstate*natoms + i] = atom->fragment[state*natoms + j];
		atom->reactive[nextstate*natoms + i] = 1;
                // New way to handle bifurcation
                // Oxygens remain available until we move on to the next hop
                toBecomeUnavailableNextHop.push_back(j);
		// Update for next cycle
		++nstates;
		++nextstate;
	      }
	    }
	  }
	  // Make donated H atom unavailable for donation now
	  atom->available[i] = 0;
	}
      }
      // Move on to next state
      ++state;
      if (hop_length > state) {
	// Set new hop level
	hop = atom->hop[state*natoms + 0];
        // Make tagged oxygens unavailable now
        for (int k=0; k<toBecomeUnavailableNextHop.size(); ++k) {
          atom->available[ toBecomeUnavailableNextHop[k] ] = 0;
        }
        // Clear the vector for next round
        toBecomeUnavailableNextHop.clear();
      }
    } // close while


    // ** Figure out which atoms belong to the environment ** //
    // Those which are environment fragments will be approximated as being equal to that
    // of the pivot state to save on computation time.
    for (int i=0; i<natoms; ++i) {
      atom->environment[i] = 1; // default to all being environment
      for (int s=0; s<nstates; ++s) {
        if (atom->reactive[s*natoms + i]) {
          atom->environment[i] = 0; 
          break;
        } 
      }
    }

    printf("State search done.\n");
    printf("Final state hop: %d\n", hop);
    printf("Number of states found: %d\n", nstates);

//#ifdef FMR_DEBUG
//#if 0
    if (fmr->print_level > 0) {
      // Debug printing
      for (int i=0; i<natoms; ++i) {
        printf("%c %12f %12f %12f", atom->symbol[i], atom->coord[3*i], atom->coord[3*i+1], atom->coord[3*i+2]); 
        for (int s=0; s<nstates; ++s) {
          printf(" %d", atom->fragment[s*natoms + i]);
        }
        printf("\n");
      }
    }
//#endif

    // Set number of states
    fmr->atom->nstates = nstates;

    // ** Check if fragment indexing has changed ** //
    for (int i=0; i<natoms*MAX_STATES; ++i) {
      if (old_fragment[i] != atom->fragment[i]) {
        fragmentIndexChange = 1; 
        break; 
      }
    }

    delete [] old_fragment;
  }

  // ** Communicate state information found by master rank ** //
  MPI_Bcast(&fmr->atom->nstates, 1, MPI_INT, MASTER_RANK, fmr->world); 
  MPI_Bcast(fmr->atom->fragment, natoms*MAX_STATES, MPI_INT, MASTER_RANK, fmr->world); 
  MPI_Bcast(fmr->atom->reactive, natoms*MAX_STATES, MPI_INT, MASTER_RANK, fmr->world); 
  MPI_Bcast(fmr->atom->environment, natoms, MPI_INT, MASTER_RANK, fmr->world); 
  MPI_Bcast(&fragmentIndexChange, 1, MPI_INT, MASTER_RANK, fmr->world); 

  // ** Check if number of states has changed ** //
  if (fmr->atom->prev_nstates != fmr->atom->nstates) {
    // Then we have to do some stuff to handle the change
    if (fmr->master_rank) {
      printf("Change in number of states. Cannot read MOs from disk.\n");
    }
    flag_state_number_change = 1;
    flag_read_MOs = 0; // Can't read MOs now b/c states appeared/disappeared
  } else {
    flag_state_number_change = 0;
    flag_read_MOs = 1; // Should be safe to read MOs (unless one state appears and one disappears 
                       // on the same step and fragment indexing changes... below)
  }
  if (fragmentIndexChange && flag_read_MOs) {
    if (fmr->master_rank) {
      printf("Change in fragment indexing. Cannot read MOs from disk.\n");
    }
    flag_read_MOs = 0; // Fragment indexing has changed since last step. Can't read MOs.
  }

}


/*-----------------------------------------------------------------
  Write Q-Chem inputs for each state's FMO calculations 
-----------------------------------------------------------------*/
void State::write_qchem_inputs(int jobtype)
{
  // Writes a separate input file for all monomers and all dimers
  // Master rank does all the work here

  if (fmr->master_rank) {

    printf("Writing Q-Chem inputs.\n");
    printf("Read MOs: %d\n", flag_read_MOs);

    Atom *atom     = fmr->atom;
    int natoms     = atom->natoms;
    int nstates    = atom->nstates;
    int nfragments = atom->nfragments;

    // ***** Loop over states ***** //
    for (int istate=0; istate<nstates; ++istate) {

      // Determine the charged reactive fragment for this state
      int chgfrag = 0;
      for (int i=0; i<natoms; ++i) {
        if (atom->reactive[istate*natoms + i]) {
          chgfrag = atom->fragment[istate*natoms + i];
          break;
        }
      }

      // Put files in directory for organization
      char state_directory[256];
      char snum[16];
      char make_directory[256];
      if (istate >= 10) {
        sprintf(snum, "%d", istate);
        sprintf(state_directory, "state_%d", istate);
      } else {
        sprintf(snum, "0%d", istate);
        sprintf(state_directory, "state_0%d", istate);
      }
      // Make the directory...
      sprintf(make_directory, "mkdir -p %s", state_directory);
      int ierr = system(make_directory);

      // *** Monomers *** //
      for (int ifrag=0; ifrag<nfragments; ++ifrag){ 
        // Get name of file to open
        char filename[256];
        if (ifrag >= 100) {
          sprintf(filename, "%s/fmo_st%sm%d.in", state_directory, snum, ifrag);
        } else if (ifrag >= 10) {
          sprintf(filename, "%s/fmo_st%sm0%d.in", state_directory, snum, ifrag);
        } else {
          sprintf(filename, "%s/fmo_st%sm00%d.in", state_directory, snum, ifrag);
        } 
        FILE *fs = fopen(filename, "w");
        if (fs == NULL) {
          char tmpstr[256];
          sprintf(tmpstr, "Failure to write Q-Chem input for file %s", filename);
          fmr->error(FLERR, tmpstr);
        }

        // Comment for labeling
        fprintf(fs, "$comment\n");
        fprintf(fs, "State %d Monomer %d\n", istate, ifrag);
        fprintf(fs, "$end\n\n");

        // $rem section
        fprintf(fs, "$rem\n");
        if (jobtype == RUN_ENERGY) 
          fprintf(fs, "jobtype sp\n");
        else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
          fprintf(fs, "jobtype force\n");
        //fprintf(fs, "exchange pbe0\n");
        //fprintf(fs, "basis 6-31+G*\n");
        fprintf(fs, "basis cc-pvdz\n");
	//fprintf(fs, "aux_basis rimp2-cc-pvdz\n");
        fprintf(fs, "exchange hf\n");
        fprintf(fs, "correlation mp2\n");
        //fprintf(fs, "correlation rimp2\n");
        fprintf(fs, "scf_convergence 6\n");
        fprintf(fs, "qm_mm true\n");
        fprintf(fs, "print_input true\n");
        fprintf(fs, "sym_ignore true\n");
        fprintf(fs, "no_reorient true\n");
        fprintf(fs, "skip_charge_self_interact 1\n");
        //fprintf(fs, "gaussian_blur true\n");
        // Read previous step MO coeffs?
        if (flag_read_MOs) fprintf(fs, "scf_guess read\n");
        fprintf(fs, "$end\n\n");

        // $molecule section
        fprintf(fs, "$molecule\n");
        if (ifrag == chgfrag) {
          fprintf(fs, "1 1\n");
        } else {
          fprintf(fs, "0 1\n");
        }
        for (int iatom=0; iatom<natoms; ++iatom) {
          if (atom->fragment[istate*natoms + iatom] == ifrag) {
            fprintf(fs, "%c %20.10lf %20.10lf %20.10lf\n",
                    atom->symbol[iatom],
                    atom->coord[3*iatom],
                    atom->coord[3*iatom+1],
                    atom->coord[3*iatom+2]
                   );
          }
        }
        fprintf(fs, "$end\n\n");

        // $external_charges section
        fprintf(fs, "$external_charges\n");
        for (int iatom=0; iatom<natoms; ++iatom) {
          if (atom->fragment[istate*natoms + iatom] != ifrag) {
            double mmq = atom->getCharge(iatom, istate);
            //fprintf(fs, "%20.10lf %20.10lf %20.10lf %16.4lf 0.12446\n",
            fprintf(fs, "%20.10lf %20.10lf %20.10lf %16.4lf\n",
                    atom->coord[3*iatom],
                    atom->coord[3*iatom+1],
                    atom->coord[3*iatom+2],
                    mmq
                   );
          }
        }
        fprintf(fs, "$end\n\n");
 
        fclose(fs);
      } // close loop over fragments for monomers


      // *** Dimers *** //
      for (int ifrag=0; ifrag<nfragments; ++ifrag){ 
        for (int jfrag=ifrag+1; jfrag<nfragments; ++jfrag){ 
	  // Get name of file to open
	  char filename[256];
          char inum[16];
          char jnum[16];
	  if (ifrag >= 100) {
            sprintf(inum, "%d", ifrag);
	  } else if (ifrag >= 10) {
            sprintf(inum, "0%d", ifrag);
	  } else {
            sprintf(inum, "00%d", ifrag);
	  } 
	  if (jfrag >= 100) {
            sprintf(jnum, "%d", jfrag);
	  } else if (jfrag >= 10) {
            sprintf(jnum, "0%d", jfrag);
	  } else {
            sprintf(jnum, "00%d", jfrag);
	  } 
	  sprintf(filename, "%s/fmo_st%sd%s-%s.in", state_directory, snum, inum, jnum);
	  FILE *fs = fopen(filename, "w");
	  if (fs == NULL) {
	    char tmpstr[256];
	    sprintf(tmpstr, "Failure to write Q-Chem input for file %s", filename);
	    fmr->error(FLERR, tmpstr);
	  }

	  // Comment for labeling
	  fprintf(fs, "$comment\n");
	  fprintf(fs, "State %d Dimer %d %d\n", istate, ifrag, jfrag);
	  fprintf(fs, "$end\n\n");

	  // $rem section
	  fprintf(fs, "$rem\n");
	  if (jobtype == RUN_ENERGY) 
	    fprintf(fs, "jobtype sp\n");
	  else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
	    fprintf(fs, "jobtype force\n");
	  //fprintf(fs, "exchange pbe0\n");
          //fprintf(fs, "basis 6-31+G*\n");
	  fprintf(fs, "basis cc-pvdz\n");
	  //fprintf(fs, "aux_basis rimp2-cc-pvdz\n");
	  fprintf(fs, "exchange hf\n");
	  fprintf(fs, "correlation mp2\n");
	  //fprintf(fs, "correlation rimp2\n");
	  fprintf(fs, "scf_convergence 6\n");
	  fprintf(fs, "qm_mm true\n");
	  fprintf(fs, "print_input true\n");
	  fprintf(fs, "sym_ignore true\n");
	  fprintf(fs, "no_reorient true\n");
          fprintf(fs, "skip_charge_self_interact 1\n");
          //fprintf(fs, "gaussian_blur true\n");
          // Read previous step MO coeffs?
          if (flag_read_MOs) fprintf(fs, "scf_guess read\n");
	  fprintf(fs, "$end\n\n");

	  // $molecule section
	  fprintf(fs, "$molecule\n");
	  if (ifrag == chgfrag || jfrag == chgfrag) {
	    fprintf(fs, "1 1\n");
	  } else {
	    fprintf(fs, "0 1\n");
	  }
	  for (int iatom=0; iatom<natoms; ++iatom) {
	    if ((atom->fragment[istate*natoms + iatom] == ifrag) ||
	        (atom->fragment[istate*natoms + iatom] == jfrag)) {
	      fprintf(fs, "%c %20.10lf %20.10lf %20.10lf\n",
		      atom->symbol[iatom],
		      atom->coord[3*iatom],
		      atom->coord[3*iatom+1],
		      atom->coord[3*iatom+2]
		     );
	    }
	  }
	  fprintf(fs, "$end\n\n");

	  // $external_charges section
	  fprintf(fs, "$external_charges\n");
	  for (int iatom=0; iatom<natoms; ++iatom) {
	    if (atom->fragment[istate*natoms + iatom] != ifrag &&
	        atom->fragment[istate*natoms + iatom] != jfrag) {
              double mmq = atom->getCharge(iatom, istate);
              //fprintf(fs, "%20.10lf %20.10lf %20.10lf %16.4lf 0.12446\n",
	      fprintf(fs, "%20.10lf %20.10lf %20.10lf %16.4lf\n",
		      atom->coord[3*iatom],
		      atom->coord[3*iatom+1],
		      atom->coord[3*iatom+2],
		      mmq
		     );
	    }
	  }
	  fprintf(fs, "$end\n\n");
   
	  fclose(fs);
        } 
      } // close loop over fragments for dimers

    } // close loop over states

    printf("Done writing Q-Chem inputs.\n");
  }

  // Hold up
  MPI_Barrier(fmr->world);
}

/*-----------------------------------------------------------------
  Update pivot state information, etc. for next step 
  Called only after matrix diagonalization such that coeffs are available
  Only the master rank needs this information since it is the only
  rank that does the state search
-----------------------------------------------------------------*/
void State::updatePivotState()
{
  int *reactive = fmr->atom->reactive;
  int *fragment = fmr->atom->fragment;
  int natoms    = fmr->atom->natoms;
  int nstates   = fmr->atom->nstates;

  if (fmr->master_rank) {
    double *GSCoeffs = fmr->matrix->GSCoeffs;

    // Next pivot state is the one with the greatest amplitude
    double max_cc = GSCoeffs[0] * GSCoeffs[0];
    int    next_pivot_state = 0;
    for (int i=1; i<nstates; ++i) { 
      double cc = GSCoeffs[i] * GSCoeffs[i];
      if (max_cc < cc) {
        max_cc = cc;
        next_pivot_state = i;
      }
    }
    printf("Pivot state for next step: %d\n", next_pivot_state);

    // Assign new reactive fragment index for pivot state
    for (int i=0; i<natoms; ++i) { 
      if (reactive[next_pivot_state*natoms + i]) {
        fmr->atom->ireactive = fragment[next_pivot_state*natoms + i];
        break;
      }
    }
    printf("Reactive fragment in next pivot state: %d\n", fmr->atom->ireactive);

    // Copy the reactive and fragment array data to state zero to serve as pivot
    // state info for next state search, but only if pivot state has changed index from being zero
    if (next_pivot_state != 0) {
      for (int i=0; i<natoms; ++i) { 
        reactive[i] = reactive[next_pivot_state*natoms + i]; 
        fragment[i] = fragment[next_pivot_state*natoms + i]; 
      }
      // Also, now it is unsafe to read MOs from disk
      flag_read_MOs = 0;      
    }     
  }

  // ** Broadcast the above to workers ** //
  MPI_Bcast(&fmr->atom->ireactive, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&next_pivot_state, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(reactive, natoms, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(fragment, natoms, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&flag_read_MOs, 1, MPI_INT, MASTER_RANK, fmr->world);
}
