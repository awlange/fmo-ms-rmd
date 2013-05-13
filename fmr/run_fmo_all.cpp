/* AWGL */
#include "fmr.h"
#include "run.h"
#include "atom.h"
#include "state.h"
#include "input.h"

#define MAX_LENGTH 1024

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Perform all the FMO calculations 
-----------------------------------------------------------------*/
void Run::do_fmo_calculations(int FORCE)
{
  // Check if we should branch to the env approximation
  if (EnvApprox) {
    do_fmo_calculations_env();
    return;
  }

  // Divide the list up across MPI ranks
  // Then, run the FMO calculations with call to serial Q-Chem in parallel

  int natoms     = fmr->atom->natoms;
  int nfragments = fmr->atom->nfragments;
  int nstates    = fmr->atom->nstates;
  int my_rank    = fmr->my_rank;
  int nf2        = nfragments*nfragments;

  n_monomers = nstates * nfragments;
  // Assuming all states have equal number of dimers, for now
  n_dimers = nstates * ((nfragments * (nfragments-1)) / 2);
  n_dimers_sq = nstates * nfragments * nfragments; // inclues self

  if (fmr->master_rank) {
    printf("Preparing to run FMO calculations:\n");
    printf("Monomer FMO calculations: %d\n", n_monomers);
    printf("Dimer FMO calculations:   %d\n", n_dimers);
  }

  // ** If number of states changed, need to deallocate memory and reallocate below ** //
  if (fmr->state->flag_state_number_change) {
    delete [] fmo_energies;
    delete [] monomer_energies;
    delete [] dimer_energies;
    fmo_energies = monomer_energies = dimer_energies = NULL;
    if (FORCE) {
      delete [] fmo_gradients;
      delete [] monomer_gradients;
      delete [] dimer_gradients;
      fmo_gradients = monomer_gradients = dimer_gradients = NULL;
    }
  }

  // ** Allocate energies and zero ** //
  if (fmo_energies == NULL)     fmo_energies     = new double [nstates];
  if (monomer_energies == NULL) monomer_energies = new double [n_monomers];
  if (dimer_energies == NULL)   dimer_energies   = new double [n_dimers_sq];
  for (int i=0; i<nstates; ++i)     fmo_energies[i]     = 0.0;
  for (int i=0; i<n_monomers; ++i)  monomer_energies[i] = 0.0;
  for (int i=0; i<n_dimers_sq; ++i) dimer_energies[i]   = 0.0;

  if (FORCE) {
    // ** Allocate gradients and zero ** //
    if (fmo_gradients == NULL)     fmo_gradients     = new double[nstates*3*natoms];
    if (monomer_gradients == NULL) monomer_gradients = new double[n_monomers*3*natoms];
    if (dimer_gradients == NULL)   dimer_gradients   = new double[n_dimers_sq*3*natoms];
    for (int i=0; i<nstates*3*natoms; ++i)     fmo_gradients[i]     = 0.0;
    for (int i=0; i<n_monomers*3*natoms; ++i)  monomer_gradients[i] = 0.0;
    for (int i=0; i<n_dimers_sq*3*natoms; ++i) dimer_gradients[i]   = 0.0;
  }

  // ** Determine load balance ** //
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
  div = n_dimers / fmr->world_size;
  rem = n_dimers % fmr->world_size;
  int ifrom_dim, ito_dim;
  if (my_rank < rem) {
    ifrom_dim = my_rank*div + my_rank;
    ito_dim   = ifrom_dim + div + 1;
  } else {  
    ifrom_dim = my_rank*div + rem;
    ito_dim   = ifrom_dim + div;
  }

  // Clock
  MPI_Barrier(fmr->world);
  double FMO_start = MPI_Wtime();
  if (fmr->master_rank) {
    printf("Running FMO calculations...\n");
  }

  // ** Make the system calls to run each FMO calculation now ** //
  char command[MAX_LENGTH];
  int ierr;
  int index_mono = 0;
  int index_dim  = 0;

  for (int istate=0; istate<nstates; ++istate) {

    char state_directory[256];
    char snum[16];
    if (istate >= 10) {
      sprintf(snum, "%d", istate);
      sprintf(state_directory, "state_%d", istate);
    } else {
      sprintf(snum, "0%d", istate);
      sprintf(state_directory, "state_0%d", istate);
    }

    // ********** Handle FMO monomers here ************ //
    for (int ifrag=0; ifrag<nfragments; ++ifrag) {
      if (ifrom_mono <= index_mono && index_mono < ito_mono) {
        char filename[256];
        char inum[16];
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
        //printf("Rank %d: %s\n", my_rank, command);

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

        if (FORCE) {
          // ** Get gradient from file ** // 
          sprintf(output_file, "%s/%s.in.gradient", state_directory, filename);
          fs = fopen(output_file, "r");
          if (fs == NULL) {
            char tmpstr[MAX_LENGTH];
            sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
            fmr->error(FLERR, tmpstr);
          }
          char line[MAX_LENGTH];
          int iatom, atnum;
          atnum = 0; // index of QM atom for storing gradient
          double gx, gy, gz;
          while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
            // Advance atnum until it matches as a QM atom index for this monomer fragment
            while ( !fmr->atom->AtomInFragment(atnum, ifrag, istate) ) {
              atnum++; 
            }
	    if ( sscanf(line, "%d %lf %lf %lf", &iatom, &gx, &gy, &gz) == 4 ) {
              monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*atnum]   = gx; 
              monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*atnum+1] = gy; 
              monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*atnum+2] = gz; 
	    }
            // Increment atnum for the next round
            atnum++;
          }
          fclose(fs);

          // ** Get field from file ** // 
          sprintf(output_file, "%s/%s.in.field", state_directory, filename);
          fs = fopen(output_file, "r");
          if (fs == NULL) {
            char tmpstr[MAX_LENGTH];
            sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
            fmr->error(FLERR, tmpstr);
          }
          atnum = 0; // index of non-QM atom for storing gradient
          while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
            // Advance atnum until it matches as a non-QM atom index for this monomer fragment
            while ( fmr->atom->AtomInFragment(atnum, ifrag, istate) ) {
              atnum++; 
            }
	    if ( sscanf(line, "%d %lf %lf %lf", &iatom, &gx, &gy, &gz) == 4 ) {
              // gx,gy,gz = the electric field
              // multiply by charge to get force (i.e. negative gradient) on atom
              double mmq = fmr->atom->getCharge(atnum, istate);
              gx *= -mmq;
              gy *= -mmq;
              gz *= -mmq;
              monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*atnum]   = gx; 
              monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*atnum+1] = gy; 
              monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*atnum+2] = gz; 
	    }
            // Increment atnum for the next round
            atnum++;
          }
          fclose(fs);
        }

      }
      ++index_mono;
    }

    // ********** Handle FMO dimers in this loop ************ //
    for (int ifrag=0; ifrag<nfragments; ++ifrag) {
      for (int jfrag=ifrag+1; jfrag<nfragments; ++jfrag) {
        if (ifrom_dim <= index_dim && index_dim < ito_dim) {
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
	  //printf("Rank %d: %s\n", my_rank, command);

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

          if (FORCE) {
            // ** Get gradient from file ** // 
            sprintf(output_file, "%s/%s.in.gradient", state_directory, filename);
            fs = fopen(output_file, "r");
            if (fs == NULL) {
              char tmpstr[MAX_LENGTH];
              sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
              fmr->error(FLERR, tmpstr);
            }
            char line[MAX_LENGTH];
            int iatom; // dummy index
            int atnum = 0; // index of QM atom for storing gradient
            double gx, gy, gz;
            while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
              // Advance atnum until it matches as a QM atom index for this dimer fragment
              while ( !(fmr->atom->AtomInFragment(atnum, ifrag, istate) || 
                        fmr->atom->AtomInFragment(atnum, jfrag, istate)) ) {
                atnum++; 
              }
  	      if ( sscanf(line, "%d %lf %lf %lf", &iatom, &gx, &gy, &gz) == 4 ) {
                // store symmetrically
                dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*atnum]   = gx; 
                dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*atnum+1] = gy; 
                dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*atnum+2] = gz; 
                dimer_gradients[(nf2*istate + nfragments*jfrag + ifrag)*3*natoms + 3*atnum]   = gx; 
                dimer_gradients[(nf2*istate + nfragments*jfrag + ifrag)*3*natoms + 3*atnum+1] = gy; 
                dimer_gradients[(nf2*istate + nfragments*jfrag + ifrag)*3*natoms + 3*atnum+2] = gz; 
	      }
              // Increment atnum for the next round
              atnum++;
            }
            fclose(fs);

            // ** Get field from file ** // 
            sprintf(output_file, "%s/%s.in.field", state_directory, filename);
            fs = fopen(output_file, "r");
            if (fs == NULL) {
              char tmpstr[MAX_LENGTH];
              sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
              fmr->error(FLERR, tmpstr);
            }
            atnum = 0; // index of non-QM atom for storing gradient
            while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
              // Advance atnum until it matches as a non-QM atom index for this dimer fragment
              while ( (fmr->atom->AtomInFragment(atnum, ifrag, istate) || 
                       fmr->atom->AtomInFragment(atnum, jfrag, istate)) ) {
                atnum++; 
              }
 	      if ( sscanf(line, "%d %lf %lf %lf", &iatom, &gx, &gy, &gz) == 4 ) {
                // gx,gy,gz = the electric field
                // multiply by charge to get force (i.e. negative gradient) on atom
                double mmq = fmr->atom->getCharge(atnum, istate);
                gx *= -mmq;
                gy *= -mmq;
                gz *= -mmq;
                // store symmetrically 
                dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*atnum]   = gx; 
                dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*atnum+1] = gy; 
                dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*atnum+2] = gz; 
                dimer_gradients[(nf2*istate + nfragments*jfrag + ifrag)*3*natoms + 3*atnum]   = gx; 
                dimer_gradients[(nf2*istate + nfragments*jfrag + ifrag)*3*natoms + 3*atnum+1] = gy; 
                dimer_gradients[(nf2*istate + nfragments*jfrag + ifrag)*3*natoms + 3*atnum+2] = gz; 
              }
              // Increment atnum for the next round
              atnum++;
            }
            fclose(fs);
          }
	}
        ++index_dim;
      }
    }

  }

  // Clock
  double FMO_end = MPI_Wtime();
  MPI_Barrier(fmr->world);
  if (fmr->master_rank) {
    printf("Finished with FMO calculations.\n");
  }

  // *** Reduce the energies from the parallel Q-Chem calls *** //
  double *rbuffer; 
  if (FORCE) rbuffer = new double [n_dimers_sq*3*natoms]; 
  else       rbuffer = new double [n_dimers_sq];

  // Monomers
  for (int i=0; i<n_monomers; ++i) rbuffer[i] = 0.0;
  MPI_Allreduce(monomer_energies, rbuffer, n_monomers, MPI_DOUBLE, MPI_SUM, fmr->world);
  for (int i=0; i<n_monomers; ++i) monomer_energies[i] = rbuffer[i];
  // Dimers
  for (int i=0; i<n_dimers_sq; ++i) rbuffer[i] = 0.0;
  MPI_Allreduce(dimer_energies, rbuffer, n_dimers_sq, MPI_DOUBLE, MPI_SUM, fmr->world);
  for (int i=0; i<n_dimers_sq; ++i) dimer_energies[i] = rbuffer[i];
  
  if (FORCE) {
    // Monomers
    for (int i=0; i<n_monomers*3*natoms; ++i) rbuffer[i] = 0.0;
    MPI_Allreduce(monomer_gradients, rbuffer, n_monomers*3*natoms, MPI_DOUBLE, MPI_SUM, fmr->world);
    for (int i=0; i<n_monomers*3*natoms; ++i) monomer_gradients[i] = rbuffer[i];
    // Dimers
    for (int i=0; i<n_dimers_sq*3*natoms; ++i) rbuffer[i] = 0.0;
    MPI_Allreduce(dimer_gradients, rbuffer, n_dimers_sq*3*natoms, MPI_DOUBLE, MPI_SUM, fmr->world);
    for (int i=0; i<n_dimers_sq*3*natoms; ++i) dimer_gradients[i] = rbuffer[i];
  }

  delete [] rbuffer;

  // *** Compute the FMO energies/forces for each state *** //
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
          double etmp = dimer_energies[nfragments*nfragments*istate + nfragments*ifrag + jfrag] - 
                        monomer_energies[nfragments*istate + ifrag] -
                        monomer_energies[nfragments*istate + jfrag];
          if (fmr->print_level > 0) {
            printf("%d--%d : %16.10f %16.10f %16.10f %16.10f\n", ifrag, jfrag, 
                    dimer_energies[nfragments*nfragments*istate + nfragments*ifrag + jfrag],
                    monomer_energies[nfragments*istate + ifrag],
                    monomer_energies[nfragments*istate + jfrag], etmp
                  );
          }
          en_fmo2 += etmp;
        }
      }
      printf("Total dimer energy:   %16.10f\n", en_fmo2);
      fmo_energies[istate] = en_fmo1 + en_fmo2;
      printf("Total FMO energy:     %16.10f\n", en_fmo1 + en_fmo2);
      printf("----------------------\n");

      // ** Do forces for this state ** //
      if (FORCE) {
        double gx, gy, gz;
        gx = gy = gz = 0.0;
        // Monomer force
        for (int ifrag=0; ifrag<nfragments; ++ifrag) {
          for (int i=0; i<natoms; ++i) {
            gx = monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i];
            gy = monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i+1];
            gz = monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i+2];
            fmo_gradients[3*natoms*istate + 3*i]   += gx;
            fmo_gradients[3*natoms*istate + 3*i+1] += gy;
            fmo_gradients[3*natoms*istate + 3*i+2] += gz;
          }
        }
        // Dimer force
        for (int ifrag=0; ifrag<nfragments; ++ifrag) {
          for (int jfrag=ifrag+1; jfrag<nfragments; ++jfrag) {
            for (int i=0; i<natoms; ++i) {
               double tmpx, tmpy, tmpz;
               gx = dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*i];  
               gy = dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*i+1];
               gz = dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*i+2];
               gx -= monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i] + 
                     monomer_gradients[(nfragments*istate + jfrag)*3*natoms + 3*i];
               gy -= monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i+1] + 
                     monomer_gradients[(nfragments*istate + jfrag)*3*natoms + 3*i+1];
               gz -= monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i+2] + 
                     monomer_gradients[(nfragments*istate + jfrag)*3*natoms + 3*i+2];
               fmo_gradients[3*natoms*istate + 3*i]   += gx;
               fmo_gradients[3*natoms*istate + 3*i+1] += gy;
               fmo_gradients[3*natoms*istate + 3*i+2] += gz;
            }
          }
        }
      }
    }
  }
  // Broadcast the FMO energy/force to worker ranks
  MPI_Bcast(fmo_energies, nstates, MPI_DOUBLE, MASTER_RANK, fmr->world);
  if (FORCE) {
    MPI_Bcast(fmo_gradients, nstates*3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
  }

#ifdef FMR_DEBUG
  if (fmr->master_rank && FORCE) {
    printf("Monomer gradients:\n");
    for (int istate=0; istate<nstates; ++istate) {
      for (int ifrag=0; ifrag<nfragments; ++ifrag) {
        printf("State %d fragment %d:\n", istate, ifrag);
        for (int i=0; i<natoms; ++i) {
          printf("%3d %12.8f %12.8f %12.8f\n", i, 
                 monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i],
                 monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i+1],
                 monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i+2]
                );
        }
      }
    }
    printf("Dimer gradients:\n");
    for (int istate=0; istate<nstates; ++istate) {
      for (int ifrag=0; ifrag<nfragments; ++ifrag) {
        for (int jfrag=ifrag+1; jfrag<nfragments; ++jfrag) {
          printf("State %d fragments %d--%d:\n", istate, ifrag, jfrag);
          for (int i=0; i<natoms; ++i) {
            printf("%3d %12.8f %12.8f %12.8f\n", i, 
                    dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*i],  
                    dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*i+1],
                    dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*i+2]
                  );
          }
        }
      }
    }
    printf("FMO gradients:\n");
    for (int istate=0; istate<nstates; ++istate) {
      printf("State %d:\n", istate);
      for (int i=0; i<natoms; ++i) {
        double gx, gy, gz;
        gx = fmo_gradients[3*natoms*istate + 3*i];   
        gy = fmo_gradients[3*natoms*istate + 3*i+1]; 
        gz = fmo_gradients[3*natoms*istate + 3*i+2];
        printf("%d : %f %f %f\n", i, gx, gy, gz);
      }
    }
  }
#endif

  if (fmr->master_rank) {
    printf(" Time for FMO: %.4f seconds\n", FMO_end - FMO_start);
    //printf(" Proc %d: Time for FMO: %.4f seconds\n", my_rank, FMO_end - FMO_start);
  }

}

