/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "input.h"
#include "run.h"
#include "dynamics.h"

#define MAX_LENGTH   256
#define SMALL_LENGTH 128

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/

Input::Input(FMR *fmr) : Pointers(fmr)
{
  read_restart = 0;
  // Default file names  
  sprintf(input_file, "%s", "fmr.in");
  sprintf(atoms_file, "%s", "fmr.atoms");
  sprintf(restart_file, "%s", "fmr.restart");
}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/

Input::~Input()
{
  
}


/*-----------------------------------------------------------------
  Read input file to get all the run specs
-----------------------------------------------------------------*/
void Input::read_input_file()
{
  // Only the master rank reads the input for FMR here
  // Specs are broadcast subsequently
  if (fmr->master_rank) {

    printf("Reading input file %s\n", input_file);

    FILE *fs = fopen(input_file, "r");
    if (fs == NULL) {
      fmr->error(FLERR, "Failure to read input file.\n");
    }

    // Go through each line of input_file until end of file
    char line[MAX_LENGTH];
    while ( fgets(line, MAX_LENGTH, fs) != NULL ) { 

      char arg0[SMALL_LENGTH];
      char arg1[SMALL_LENGTH];

      if ( sscanf(line, "%s %s", arg0, arg1) == 2 ) {
        // Two argument line
     
	if ( strcmp(arg0, "QChem") == 0 ) {
	  printf("QChem executable: %s\n", arg1);
	  sprintf(fmr->run->qchem_exec,"%s", arg1); 
	}
	else if ( strcmp(arg0, "Scratch") == 0 ) {
	  printf("QChem scratch: %s\n", arg1);
	  sprintf(fmr->run->qchem_scratch,"%s", arg1); 
	}
	else if ( strcmp(arg0, "Atoms_File") == 0 ) {
	  printf("Atoms file: %s\n", arg1);
	  sprintf(atoms_file,"%s", arg1); 
	}
	else if ( strcmp(arg0, "Restart_File") == 0 ) {
	  printf("Restart file: %s\n", arg1);
	  sprintf(restart_file,"%s", arg1); 
	}
	else if ( strcmp(arg0, "Run_Type") == 0 ) {
          printf("Run type: %s\n", arg1);
	  if      ( strcmp(arg1, "energy") == 0) fmr->run->run_type = RUN_ENERGY;
	  else if ( strcmp(arg1, "force" ) == 0) fmr->run->run_type = RUN_FORCE;
	  else if ( strcmp(arg1, "moldyn") == 0) fmr->run->run_type = RUN_MOLDYN;

	  if (fmr->run->run_type < RUN_TYPE_MIN || fmr->run->run_type > RUN_TYPE_MAX)
	    fmr->error(FLERR, "Selected run type is not understood or does not exist.");
        }
	else if ( strcmp(arg0, "Max_Hops") == 0 ) {
          fmr->state->max_hops = atoi(arg1);
	  printf("Max hops: %d\n", fmr->state->max_hops);
	}
	else if ( strcmp(arg0, "PrintLevel") == 0 ) {
          fmr->print_level = atoi(arg1);
	}
        // Molecular dynamics options
        else if ( strcmp(arg0, "nTimeSteps") == 0) {
          fmr->dynamics->nTimeSteps = atoi(arg1); 
        }
        else if ( strcmp(arg0, "timestep") == 0) {
          fmr->dynamics->dt = atof(arg1); // Read in as femtoseconds
          fmr->dynamics->dt *= fmr->math->fem2au; // Convert to atomic units of time 
        }
        else if ( strcmp(arg0, "vOption") == 0) {
	  if      ( strcmp(arg1, "zero")  == 0) fmr->dynamics->vOption = DYNAMICS_ZERO;
	  else if ( strcmp(arg1, "read")  == 0) fmr->dynamics->vOption = DYNAMICS_READ;
	  else if ( strcmp(arg1, "rand")  == 0) fmr->dynamics->vOption = DYNAMICS_RAND;

	  if (fmr->dynamics->vOption < DYNAMICS_MIN || fmr->dynamics->vOption > DYNAMICS_MAX)
	    fmr->error(FLERR, "Selected vOption is not understood or does not exist.");
        }
        else if ( strcmp(arg0, "Thermostat") == 0) {
	  if      ( strcmp(arg1, "off")       == 0) fmr->dynamics->iThermostat = THERMOSTAT_OFF;
	  else if ( strcmp(arg1, "langevin")  == 0) fmr->dynamics->iThermostat = THERMOSTAT_LANGEVIN;
	  else if ( strcmp(arg1, "berendsen") == 0) fmr->dynamics->iThermostat = THERMOSTAT_BERENDSEN;
          else    fmr->error(FLERR, "Thermostat option unrecognized.");
        }
        else if ( strcmp(arg0, "tau") == 0) {
          fmr->dynamics->tau = atof(arg1); // in femtoseconds
          fmr->dynamics->tau *= fmr->math->fem2au; // Convert to a.u.
        }
        else if ( strcmp(arg0, "Temperature") == 0) {
          fmr->dynamics->targetTemperature = atof(arg1); // in Kelvin 
        }
        else if ( strcmp(arg0, "vseed") == 0) {
          fmr->dynamics->vseed = atoi(arg1); 
        }
	else if ( strcmp(arg0, "Traj_File") == 0 ) {
	  printf("Trajectory file: %s\n", arg1);
	  sprintf(fmr->dynamics->trajFile,"%s", arg1); 
	}
        else if ( strcmp(arg0, "BB") == 0) {
          printf("Using BB param: %f\n", atof(arg1));
          fmr->math->BB = atof(arg1); 
        }

      }
      else if ( sscanf(line, "%s", arg0) == 1 ) {
        // One argument line

	if ( strcmp(arg0, "FMO_only") == 0 ) {
          printf("FMO only calcuation detected.\n");
	  fmr->run->FMO_only = 1;
        }
	else if ( strcmp(arg0, "EnvApprox") == 0 ) {
	  printf("Environment approximation turned on.\n");
          fmr->run->EnvApprox = 1;
	}
	else if ( strcmp(arg0, "Restart") == 0 ) {
	  printf("Will try to restart run from file.\n");
          read_restart = 1;
          // Set the appropriate velocity option
          fmr->dynamics->vOption = DYNAMICS_READ;
	}

      }
    }

    fclose(fs);

  }
  MPI_Barrier(fmr->world);

  // Communicate run specs to all other ranks
  MPI_Bcast(&fmr->run->run_type, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->run->qchem_exec, MAX_LENGTH, MPI_CHAR, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->run->qchem_scratch, MAX_LENGTH, MPI_CHAR, MASTER_RANK, fmr->world);
  MPI_Bcast(&atoms_file, MAX_LENGTH, MPI_CHAR, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->run->FMO_only, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->run->EnvApprox, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->state->max_hops, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->print_level, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->dynamics->nTimeSteps, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->dynamics->dt, 1, MPI_DOUBLE, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->dynamics->vOption, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->dynamics->targetTemperature, 1, MPI_DOUBLE, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->dynamics->vseed, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->dynamics->trajFile, MAX_LENGTH, MPI_CHAR, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->dynamics->iThermostat, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->dynamics->tau, 1, MPI_DOUBLE, MASTER_RANK, fmr->world);
  MPI_Bcast(&read_restart, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->math->BB, 1, MPI_DOUBLE, MASTER_RANK, fmr->world);

}


/*-----------------------------------------------------------------
  Read atoms file to get initial atom info
-----------------------------------------------------------------*/
void Input::read_atoms_file()
{
  // Only the master rank reads the atoms file 
  // Info is broadcast subsequently
  if (fmr->master_rank) {

    printf("Reading atoms file %s\n", atoms_file);

    FILE *fs = fopen(atoms_file, "r");
    if (fs == NULL) {
      fmr->error(FLERR, "Failure to read atoms file.\n");
    }

    char line[MAX_LENGTH];
    int line_index = 0;
    int atom_index = 0;
    // Go through each line of until end of file
    while ( fgets(line, MAX_LENGTH, fs) != NULL ) { 

      if (line_index == 0) {
        // First line syntax: 
        // (# atoms) (# fragments) (index of reactive fragment)
        sscanf(line, "%d %d %d", &fmr->atom->natoms, 
	       &fmr->atom->nfragments, &fmr->atom->ireactive);

        int natoms = fmr->atom->natoms;
	// Based on this data, allocate atom arrays
        fmr->atom->coord = new double [3*natoms];
        fmr->atom->force = new double [3*natoms];
        fmr->atom->veloc = new double [3*natoms];
        fmr->atom->mass  = new double [natoms];
	fmr->atom->symbol    = new char [natoms];
	fmr->atom->fragment  = new int [natoms*MAX_STATES];
	fmr->atom->reactive  = new int [natoms*MAX_STATES];
	fmr->atom->available = new int [natoms];
	fmr->atom->hop       = new int [natoms*MAX_STATES];
        fmr->atom->environment  = new int [natoms];
        // Initialize fragment array
        for (int k=0; k<natoms*MAX_STATES; ++k) {
          fmr->atom->fragment[k] = -1;
        }

      } else {
        // Error check
        if (atom_index > fmr->atom->natoms && atom_index > 0) {
          sprintf(line, "Error reading atoms. natoms = %d atom_index = %d", 
                  fmr->atom->natoms, atom_index);
          fmr->error(FLERR, line); 
        }
        // Regular line syntax:
	// (Atom symbol) (X) (Y) (Z) (Fragment index in state 0, pivot state)
        // Coordinates are in Angstroms
        char tmpstr[MAX_LENGTH];
        sscanf(line, "%s %lf %lf %lf %d", 
               tmpstr,  
	       &(fmr->atom->coord[3*atom_index]),
	       &(fmr->atom->coord[3*atom_index+1]),
	       &(fmr->atom->coord[3*atom_index+2]),
	       &(fmr->atom->fragment[atom_index])
	      );
        // remove whitespace from tmpstr to get the atom symbol character
        for (int i=0; i<MAX_LENGTH; ++i) {
          if (tmpstr[i] != ' ') {
            fmr->atom->symbol[atom_index] = tmpstr[i];
            break;
          }
        }
	++atom_index;
      }
      ++line_index;
    }

    fclose(fs);
  }
  MPI_Barrier(fmr->world);


  // ** Communicate data to other ranks ** //
  MPI_Bcast(&fmr->atom->natoms, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->atom->nfragments, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->atom->ireactive, 1, MPI_INT, MASTER_RANK, fmr->world);
  int natoms = fmr->atom->natoms;
  if (!fmr->master_rank) {
    // Based on broadcast data, allocate atom arrays
    fmr->atom->coord = new double [3*natoms];
    fmr->atom->force = new double [3*natoms];
    fmr->atom->veloc = new double [3*natoms];
    fmr->atom->mass  = new double [natoms];
    fmr->atom->symbol    = new char [natoms];
    fmr->atom->fragment  = new int [natoms*MAX_STATES];
    fmr->atom->reactive  = new int [natoms*MAX_STATES];
    fmr->atom->environment  = new int [natoms];
    // I don't think the worker ranks need this stuff, so save memory
    // It's specific to the state search, which only master rank does
    //fmr->atom->available = new int [natoms];
    //fmr->atom->hop       = new int [natoms*MAX_STATES];

    // Initialize fragment array
    for (int k=0; k<natoms*MAX_STATES; ++k) {
      fmr->atom->fragment[k] = -1;
    }
  }
  MPI_Bcast(fmr->atom->symbol, natoms, MPI_CHAR, MASTER_RANK, fmr->world);
  MPI_Bcast(fmr->atom->coord, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
  MPI_Bcast(fmr->atom->fragment, natoms, MPI_INT, MASTER_RANK, fmr->world);

  // Set the atomic masses
  fmr->atom->setAtomMasses();

  if(fmr->master_rank) {
    printf("Read in %d atoms with %d fragments.\n", natoms, fmr->atom->nfragments);
  }

}

/*-----------------------------------------------------------------
  Write restart file. Always overwrites.
  Syntax:
  [Number of atoms] [Number of fragments] [Reactive fragment in pivot state] (1st line only)
  StepNumber (2nd line only) 
  [Atomic Symbol] X Y Z [Pivot state fragment IDs] VX VY VZ
-----------------------------------------------------------------*/
void Input::write_restart_file()
{
  if (fmr->master_rank) {
    int natoms    = fmr->atom->natoms;
    double *coord = fmr->atom->coord;
    double *veloc = fmr->atom->veloc;
    FILE *fs = fopen(restart_file, "w");
    if (fs == NULL) {
      char tmpstr[256];
      sprintf(tmpstr, "Failure to write to restart file %s", restart_file);
      fmr->error(FLERR, tmpstr);
    }
    fprintf(fs, "%d %d %d\n", natoms, fmr->atom->nfragments, fmr->atom->ireactive);
    fprintf(fs, "%d\n", fmr->dynamics->iCurrentStep);
    for (int i=0; i<natoms; ++i) {
      fprintf(fs, "%c %16.12f %16.12f %16.12f %3d %16.12f %16.12f %16.12f\n", 
                  fmr->atom->symbol[i],
                  coord[3*i], coord[3*i+1], coord[3*i+2],
                  fmr->atom->fragment[i], // assuming next pivot state is state index 0, as in updatePivotState
                  veloc[3*i], veloc[3*i+1], veloc[3*i+2]
             );
    }
    fclose(fs);
  }
}

/*-----------------------------------------------------------------
  Read restart file instead of reading atoms file. 
  Syntax:
  [Number of atoms] [Number of fragments] [Reactive fragment in pivot state] (1st line only)
  StepNumber (2nd line only) 
  [Atomic Symbol] X Y Z [Pivot state fragment IDs] VX VY VZ
-----------------------------------------------------------------*/
void Input::read_restart_file()
{
  if (fmr->master_rank) {
    printf("Reading restart data from file %s\n", restart_file);    

    FILE *fs = fopen(restart_file, "r");
    if (fs == NULL) {
      char tmpstr[256];
      sprintf(tmpstr, "Failure to read from restart file %s", restart_file);
      fmr->error(FLERR, tmpstr);
    }

    char line[MAX_LENGTH];
    int line_index = 0;
    int atom_index = 0;
    // Go through each line until end of file
    while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
      if (line_index == 0) {
        // First line 
        sscanf(line, "%d %d %d", &fmr->atom->natoms, &fmr->atom->nfragments, &fmr->atom->ireactive);
	// Based on this data, allocate atom arrays
        int natoms = fmr->atom->natoms;
        fmr->atom->coord = new double [3*natoms];
        fmr->atom->force = new double [3*natoms];
        fmr->atom->veloc = new double [3*natoms];
        fmr->atom->mass  = new double [natoms];
	fmr->atom->symbol    = new char [natoms];
	fmr->atom->fragment  = new int [natoms*MAX_STATES];
	fmr->atom->reactive  = new int [natoms*MAX_STATES];
	fmr->atom->available = new int [natoms];
	fmr->atom->hop       = new int [natoms*MAX_STATES];
        fmr->atom->environment  = new int [natoms];
        // Initialize fragment array
        for (int k=0; k<natoms*MAX_STATES; ++k) {
          fmr->atom->fragment[k] = -1;
        }
      } 
      else if (line_index == 1) {
        // Second line 
        sscanf(line, "%d", &fmr->dynamics->iCurrentStep);
      }
      else { 
        int i = atom_index;
        double *coord = fmr->atom->coord;
        double *veloc = fmr->atom->veloc;
        sscanf(line, "%c %lf %lf %lf %d %lf %lf %lf\n", 
                  &fmr->atom->symbol[i],
                  &coord[3*i], &coord[3*i+1], &coord[3*i+2],
                  &fmr->atom->fragment[i], // assuming next pivot state is state index 0, as in updatePivotState
                  &veloc[3*i], &veloc[3*i+1], &veloc[3*i+2]
              );
        atom_index++;
      }
      line_index++;
    }
    fclose(fs);
  }
  MPI_Barrier(fmr->world);

  // ** Communicate data to other ranks ** //
  MPI_Bcast(&fmr->dynamics->iCurrentStep, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->atom->natoms, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->atom->nfragments, 1, MPI_INT, MASTER_RANK, fmr->world);
  MPI_Bcast(&fmr->atom->ireactive, 1, MPI_INT, MASTER_RANK, fmr->world);
  int natoms = fmr->atom->natoms;
  if (!fmr->master_rank) {
    // Based on broadcast data, allocate atom arrays
    fmr->atom->coord = new double [3*natoms];
    fmr->atom->force = new double [3*natoms];
    fmr->atom->veloc = new double [3*natoms];
    fmr->atom->mass  = new double [natoms];
    fmr->atom->symbol    = new char [natoms];
    fmr->atom->fragment  = new int [natoms*MAX_STATES];
    fmr->atom->reactive  = new int [natoms*MAX_STATES];
    fmr->atom->environment  = new int [natoms];
    // I don't think the worker ranks need this stuff, so save memory
    // It's specific to the state search, which only master rank does
    //fmr->atom->available = new int [natoms];
    //fmr->atom->hop       = new int [natoms*MAX_STATES];

    // Initialize fragment array
    for (int k=0; k<natoms*MAX_STATES; ++k) {
      fmr->atom->fragment[k] = -1;
    }
  }
  MPI_Bcast(fmr->atom->symbol, natoms, MPI_CHAR, MASTER_RANK, fmr->world);
  MPI_Bcast(fmr->atom->coord, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
  MPI_Bcast(fmr->atom->veloc, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
  MPI_Bcast(fmr->atom->fragment, natoms, MPI_INT, MASTER_RANK, fmr->world);


  // Set the atomic masses
  fmr->atom->setAtomMasses();

  if(fmr->master_rank) {
    printf("Read in %d atoms with %d fragments.\n", natoms, fmr->atom->nfragments);
  }
}
