/* AWGL */
#include "fmr.h"
#include "run.h"
#include "atom.h"
#include "state.h"
#include "input.h"

#define MAX_LENGTH 1024

using namespace FMR_NS;


/*-----------------------------------------------------------------
  Perform FMO calculations but use the env approximation to
  cut down on dimer similar calculations 
-----------------------------------------------------------------*/
void Run::do_fmo_calculations_env()
{

  if (fmr->master_rank) printf("Applying the environment dimer approximation.\n");

  // Divide the list up across MPI ranks
  // Then, run the FMO calculations with call to serial Q-Chem in parallel

  int natoms     = fmr->atom->natoms;
  int nfragments = fmr->atom->nfragments;
  int nstates    = fmr->atom->nstates;
  int my_rank    = fmr->my_rank;

  // Miscellaneous
  char command[MAX_LENGTH];
  int ierr;
  int index_mono = 0;
  int index_dim  = 0;
  int istate;
  char state_directory[256];
  char snum[16];
  char filename[256];
  char inum[16];
  char jnum[16];

  // Number of fragment calculations
  n_monomers = nstates * nfragments;
  // Assuming all states have equal number of dimers, for now
  n_dimers = nstates * ((nfragments * (nfragments-1)) / 2);
  n_dimers_sq = nstates * nfragments * nfragments; // includes self

  if (fmr->master_rank) {
    printf("Preparing to run FMO calculations:\n");
    printf("Monomer FMO calculations: %d\n", n_monomers);
    printf("Dimer FMO calculations:   %d\n", n_dimers);
  }

  // ** Make list of environment fragments ** //
  int sum_frag_env = 0;
  int *frag_env_list = new int [nfragments];
  for (int ifrag=0; ifrag<nfragments; ++ifrag) {
    frag_env_list[ifrag] = 1;
    for (int i=0; i<natoms && frag_env_list[ifrag]; ++i) {
      if (!atom->environment[i]) {
	for (int s=0; s<nstates; ++s) {
	  if (atom->fragment[s*natoms + i] == ifrag) {
	    frag_env_list[ifrag] = 0;
	    break;
	  }
	}
      }
    }
    sum_frag_env += frag_env_list[ifrag];
  }

  // ** Adjust number of dimers based on the env dimers ** //
  int n_env_dim = (nstates-1)*sum_frag_env*(sum_frag_env-1)/2;
  n_dimers -= n_env_dim;
  if (fmr->master_rank) {
    printf("Environment list:\n");
    for (int ifrag=0; ifrag<nfragments; ++ifrag) {
      printf("Frag %3d: %d\n", ifrag, frag_env_list[ifrag]);
    }
    printf("Env fragments: %d\n", sum_frag_env);
    printf("Env dimers:    %d\n", n_env_dim); 
    printf("Unique dimers: %d\n", n_dimers);
  }

  // Clock
  MPI_Barrier(fmr->world);
  double FMO_start = MPI_Wtime();
  if (fmr->master_rank) {
    printf("Running FMO calculations...\n");
  }


  // ** Allocate energies and zero ** //
  if (fmo_energies == NULL)
    fmo_energies = new double [nstates];
  if (monomer_energies == NULL)
    monomer_energies = new double [n_monomers];
  if (dimer_energies == NULL)
    dimer_energies   = new double [n_dimers_sq];
  for (int i=0; i<nstates; ++i)     fmo_energies[i]     = 0.0;
  for (int i=0; i<n_monomers; ++i)  monomer_energies[i] = 0.0;
  for (int i=0; i<n_dimers_sq; ++i) dimer_energies[i]   = 0.0;
  double *rbuffer = new double [n_dimers_sq];

  // ** Determine load balance for all monomers ** //
  int div = n_monomers / fmr->world_size;
  int rem = n_monomers % fmr->world_size;
  int ifrom_mono, ito_mono;
  if (my_rank < rem) {
    ifrom_mono = my_rank*div + my_rank;
    ito_mono   = ifrom_mono + div + 1;
  } else {  
    ifrom_mono = my_rank*div + rem;
    ito_mono   = ifrom_mono + div;
  }


  // ********** Handle all state FMO monomers here ************ //
  for (istate=0; istate<nstates; ++istate) {
    if (istate >= 10) {
      sprintf(snum, "%d", istate);
      sprintf(state_directory, "state_%d", istate);
    } else {
      sprintf(snum, "0%d", istate);
      sprintf(state_directory, "state_0%d", istate);
    }
    for (int ifrag=0; ifrag<nfragments; ++ifrag) {
      if (ifrom_mono <= index_mono && index_mono < ito_mono) {
        if (ifrag >= 100) {
          sprintf(inum, "%d", ifrag);
        } else if (ifrag >= 10) {
          sprintf(inum, "0%d", ifrag);
        } else {
          sprintf(inum, "00%d", ifrag);
        }
        sprintf(filename, "fmo_st%sm%s", snum, inum);
        sprintf(command, "%s %s/%s.in %s/%s/ > %s/%s.out", 
                qchem_exec,
                state_directory,
                filename,
                qchem_scratch,
                filename,
                state_directory,
                filename
               );
        // ** The system call ** //
        ierr = system(command);
        // ** Check for error ** //
        if (ierr) {
          printf("Q-Chem run error on rank %d:\n", fmr->my_rank);
          fmr->error(FLERR, command);
        }
        // ** Open output file and get the energy ** //
        char output_file[MAX_LENGTH];
        sprintf(output_file, "%s/%s.in.energy", state_directory, filename);
        FILE *fs = fopen(output_file, "r");
        if (fs == NULL) {
          char tmpstr[MAX_LENGTH];
          sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
          fmr->error(FLERR, tmpstr);
        }
        char line[MAX_LENGTH];
        double en;
        while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
	  if ( sscanf(line, "%lf", &en) == 1 ) {
            monomer_energies[nfragments*istate + ifrag] = en; 
	  }
        }
        fclose(fs);
      }
      ++index_mono;
    }
  }
  
  //printf("Proc %d: Done with my monomers.\n", my_rank);


  // ** Determine load balance for pivot state dimers ** //
  int nn = nfragments * (nfragments-1) / 2;
  div = nn / fmr->world_size;
  rem = nn % fmr->world_size;
  int ifrom_dim, ito_dim;
  if (my_rank < rem) {
    ifrom_dim = my_rank*div + my_rank;
    ito_dim   = ifrom_dim + div + 1;
  } else {  
    ifrom_dim = my_rank*div + rem;
    ito_dim   = ifrom_dim + div;
  }

  // ********** Handle pivot state FMO dimers in this loop ************ //
  istate = 0;
  sprintf(snum, "0%d", istate);
  sprintf(state_directory, "state_0%d", istate);


    for (int ifrag=0; ifrag<nfragments; ++ifrag) {
      for (int jfrag=ifrag+1; jfrag<nfragments; ++jfrag) {
        if (ifrom_dim <= index_dim && index_dim < ito_dim) {
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
          sprintf(filename, "fmo_st%sd%s-%s", snum, inum, jnum);
	  sprintf(command, "%s %s/%s.in %s/%s/ > %s/%s.out", 
		  qchem_exec,
                  state_directory,
		  filename,
		  qchem_scratch,
		  filename,
                  state_directory,
		  filename
		 );
          // ** The system call ** //
	  ierr = system(command);
          // ** Check for error ** //
          if (ierr) {
            printf("Q-Chem run error on rank %d:\n", fmr->my_rank);
            fmr->error(FLERR, command);
          }
	  // ** Open output file and get the energy ** //
          char output_file[MAX_LENGTH];
          sprintf(output_file, "%s/%s.in.energy", state_directory, filename);
	  FILE *fs = fopen(output_file, "r");
	  if (fs == NULL) {
	    char tmpstr[MAX_LENGTH];
	    sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
	    fmr->error(FLERR, tmpstr);
	  }
	  char line[MAX_LENGTH];
	  double en;
	  while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
	    if ( sscanf(line, "%lf", &en) == 1 ) {
              // save symmetrized
	      dimer_energies[nfragments*ifrag + jfrag] = en; 
	      dimer_energies[nfragments*jfrag + ifrag] = en; 
	    }
	  }
	  fclose(fs);
	}
        ++index_dim;
      }
    }

  //printf("Proc %d: Done with my pivot dimers.\n", my_rank);


  // Unnecessary to synchronize here. Keep going!


  // *** Redetermine load balance for remaining states *** //
  nn = 0;
  for (int ifrag=0; ifrag<nfragments; ++ifrag) {
    for (int jfrag=ifrag+1; jfrag<nfragments; ++jfrag) {
      if ( !(frag_env_list[ifrag] && frag_env_list[jfrag]) ) nn++;
    }
  }
  nn *= (nstates-1);
  div = nn / fmr->world_size;
  rem = nn % fmr->world_size;
  if (my_rank < rem) {
    ifrom_dim = my_rank*div + my_rank;
    ito_dim   = ifrom_dim + div + 1;
  } else {  
    ifrom_dim = my_rank*div + rem;
    ito_dim   = ifrom_dim + div;
  }
  index_dim = 0;
  //printf("Proc %d: ifrom = %d ito = %d tot = %d.\n", my_rank, ifrom_dim, ito_dim, ito_dim-ifrom_dim);

  // ********** Handle rest of FMO dimers in this loop, approximating the env-env ************ //
  for (istate=1; istate<nstates; ++istate) {

    //printf("Proc %d: On state %d of other dimers.\n", my_rank, istate);

    if (istate >= 10) {
      sprintf(snum, "%d", istate);
      sprintf(state_directory, "state_%d", istate);
    } else {
      sprintf(snum, "0%d", istate);
      sprintf(state_directory, "state_0%d", istate);
    }

    for (int ifrag=0; ifrag<nfragments; ++ifrag) {
      for (int jfrag=ifrag+1; jfrag<nfragments; ++jfrag) {
        if (frag_env_list[ifrag] && frag_env_list[jfrag]) {
          // Do nothing! Leave as zero. Reduction below.
        } else {
	  if (ifrom_dim <= index_dim && index_dim < ito_dim) {
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
	    sprintf(filename, "fmo_st%sd%s-%s", snum, inum, jnum);
	    sprintf(command, "%s %s/%s.in %s/%s/ > %s/%s.out", 
		    qchem_exec,
		    state_directory,
		    filename,
		    qchem_scratch,
		    filename,
		    state_directory,
		    filename
		   );
	    // ** The system call ** //
	    ierr = system(command);
	    // ** Check for error ** //
	    if (ierr) {
	      printf("Q-Chem run error on rank %d:\n", fmr->my_rank);
	      fmr->error(FLERR, command);
	    }
	    // ** Open output file and get the energy ** //
	    char output_file[MAX_LENGTH];
	    //sprintf(output_file, "%s/%s.out", state_directory, filename);
	    sprintf(output_file, "%s/%s.in.energy", state_directory, filename);
	    FILE *fs = fopen(output_file, "r");
	    if (fs == NULL) {
	      char tmpstr[MAX_LENGTH];
	      sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
	      fmr->error(FLERR, tmpstr);
	    }
	    char line[MAX_LENGTH];
	    double en;
	    while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
	      if ( sscanf(line, "%lf", &en) == 1 ) {
		// save symmetrized
		dimer_energies[nfragments*nfragments*istate + nfragments*ifrag + jfrag] = en; 
		dimer_energies[nfragments*nfragments*istate + nfragments*jfrag + ifrag] = en; 
	      }
	    }
	    fclose(fs);
	  }
          // Index non-env dimers here
          ++index_dim;
        }
      }
    }

  }
  //printf("Proc %d: Done with my other dimers.\n", my_rank);

  // Clock
  double FMO_end = MPI_Wtime();
  MPI_Barrier(fmr->world);
  if (fmr->master_rank) {
    printf("Finished with FMO calculations.\n");
  }

  // *** Reduce the energies from the parallel Q-Chem calls *** //
  // -- Monomers
  for (int i=0; i<n_monomers; ++i) rbuffer[i] = 0.0;
  MPI_Reduce(monomer_energies, rbuffer, n_monomers, MPI_DOUBLE, MPI_SUM, MASTER_RANK, fmr->world);
  if (fmr->master_rank) for (int i=0; i<n_monomers; ++i) monomer_energies[i] = rbuffer[i];
  // -- Dimers 
  nn = n_dimers_sq; 
  for (int i=0; i<nn; ++i) rbuffer[i] = 0.0;
  MPI_Reduce(dimer_energies, rbuffer, nn, MPI_DOUBLE, MPI_SUM, MASTER_RANK, fmr->world);
  if (fmr->master_rank) for (int i=0; i<nn; ++i) dimer_energies[i] = rbuffer[i]; 

  // *** Compute the FMO energies for each state *** //
  if (fmr->master_rank) {
    for (int istate=0; istate<nstates; ++istate) {
      printf("----- State %4d -----\n", istate);
      if (fmr->print_level > 0) {
        printf("Monomer : EI\n");
      }
      double en_fmo1 = 0.0;
      for (int ifrag=0; ifrag<nfragments; ++ifrag) {
        if (fmr->print_level > 0) {
          printf("%4d : %16.10f\n", ifrag, monomer_energies[nfragments*istate + ifrag]);
        }
        en_fmo1 += monomer_energies[nfragments*istate + ifrag];
      }
      printf("Total monomer energy: %16.10f\n", en_fmo1);
      if (fmr->print_level > 0) {
        printf("Dimer : EIJ EI EJ (EIJ - EI - EJ)\n");
      }
      double en_fmo2 = 0.0;
      for (int ifrag=0; ifrag<nfragments; ++ifrag) {
        for (int jfrag=ifrag+1; jfrag<nfragments; ++jfrag) {
          double edim  = dimer_energies[nfragments*nfragments*istate + nfragments*ifrag + jfrag];
          double emoni = monomer_energies[nfragments*istate + ifrag];
          double emonj = monomer_energies[nfragments*istate + jfrag];
          if (frag_env_list[ifrag] && frag_env_list[jfrag]) {
            // Use pivot state energy instead!
            edim  = dimer_energies[nfragments*ifrag + jfrag]; 
            emoni = monomer_energies[ifrag];
            emonj = monomer_energies[jfrag];
          }
          double etmp = edim - emoni - emonj; 
          if (fmr->print_level > 0) {
            printf("%d--%d : %16.10f %16.10f %16.10f %16.10f\n", ifrag, jfrag, 
                    edim, emoni, emonj, etmp); 
          }
          en_fmo2 += etmp;
        }
      }
      printf("Total dimer energy:   %16.10f\n", en_fmo2);
      fmo_energies[istate] = en_fmo1 + en_fmo2;
      printf("Total FMO energy:     %16.10f\n", en_fmo1 + en_fmo2);
      printf("----------------------\n");
    }
  }
  // Broadcast the FMO energy to worker ranks
  MPI_Bcast(fmo_energies, nstates, MPI_DOUBLE, MASTER_RANK, fmr->world);

  if (fmr->master_rank) {
    printf(" Time for FMO: %.4f seconds\n", FMO_end - FMO_start);
    //printf(" Proc %d: Time for FMO: %.4f seconds\n", my_rank, FMO_end - FMO_start);
  }

  // Memory clean up
  delete [] rbuffer;
  delete [] frag_env_list;
}

