/* AWGL */
#ifndef FMR_H
#define FMR_H

// NOTE: FMR is short for FMO-MS-RMD, which is short for:
// Fragment Molecular Orbital Multi-state Reactive Molecular Dynamics

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

namespace FMR_NS {

#define FLERR __FILE__,__LINE__
#define MASTER_RANK 0

class FMR {
 public:

  // ** Constructor/destructor ** //
  FMR(int, char **, MPI_Comm);  
  ~FMR();                       

  // ** MPI business ** //
  MPI_Comm world;               // MPI communicator for world
  int world_size;               // Size of world
  int my_rank;                  // My rank in world

  // ** Class pointers ** //
  class Math *math;             // Handles math stuff, etc. 
  class Atom *atom;             // Handles atom information 
  class Input *input;           // Handles input file and atoms file
  class Run *run;               // Handles run specs and manages calculations
  class Matrix *matrix;         // Handles Hamiltonian matrix, etc.
  class State *state;           // Handles state search algorithm and QChem input generation
  class Dynamics *dynamics;     // Handles molecular dynamics

  // ** Variables ** //
  int master_rank;
  int print_level;

  // ** Functions ** //
  void parse_command_line(int, char **);
  void execute();                              // Executes the calculation
  void error(const char *, int, const char *); // Controlled error crash

};

}

#endif
