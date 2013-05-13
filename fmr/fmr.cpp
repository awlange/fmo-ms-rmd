/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "input.h"
#include "run.h"
#include "state.h"
#include "matrix.h"
#include "dynamics.h"

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/

FMR::FMR(int argc, char **argv, MPI_Comm communicator)
{

  // Handle MPI stuff
  world = communicator;
  MPI_Comm_size(world, &world_size);
  MPI_Comm_rank(world, &my_rank);
  master_rank = (my_rank == MASTER_RANK) ? 1 : 0;
  if (master_rank) {
    printf("\n----------------------------------------------\n");
    printf("               FMO-MS-RMD\n");
    printf("          Author: Adrian W. Lange\n"); 
    printf("----------------------------------------------\n");
    printf("FMO-MS-RMD on %d MPI ranks.\n\n", world_size);
  }

  // Create some of the basic pointers
  math     = new Math(this);
  atom     = new Atom(this);
  input    = new Input(this);
  run      = new Run(this);
  matrix   = new Matrix(this);
  state    = new State(this);
  dynamics = new Dynamics(this);

  print_level = 0;
  
  // Parse the command line input
  parse_command_line(argc, argv);

}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/

FMR::~FMR()
{
  if (matrix   != NULL) delete matrix;
  if (input    != NULL) delete input; 
  if (run      != NULL) delete run; 
  if (state    != NULL) delete state; 
  if (dynamics != NULL) delete dynamics;
  // NOTE: atom destrcutor must come last b/c other destructors need to know it's info in their destructors 
  if (atom     != NULL) delete atom; 
}


/*-----------------------------------------------------------------
  Parse command line input 
-----------------------------------------------------------------*/
void FMR::parse_command_line(int argc, char **argv)
{
  int print_help = 0;
  if (argc < 1) {
    print_help = 1;
  }
  else {
    if (argv[1] == "-h" || argv[1] == "h") {
      print_help = 1;
    }
  }

  // Crash with help menu if detected error
  if (print_help) {
    if (master_rank) {
      printf("----- FMR Help -----\n");
      printf("Command syntax:\n");
      printf("fmr.exe [input filename]\n");
      printf("\n\n");
      printf("Exiting...\n");
    }
    MPI_Barrier(world);
    MPI_Finalize();
    exit(1);
  }

  if (argc > 1) {
    // Get input file name and store it 
    sprintf(input->input_file, "%s", argv[1]);
  }
  else {
    // Assume default file names
    if (master_rank) {
      printf("Assuming default file names.\n");
    }
  }

}

/*-----------------------------------------------------------------
  Crash with error message 
-----------------------------------------------------------------*/

void FMR::error(const char *file, int line, const char *str)
{
  // Controlled error crash
  printf("Error on rank: %d   File: %s   Line: %d   Message: %s\n", 
         my_rank, file, line, str);
  MPI_Abort(world,1);
}

/*-----------------------------------------------------------------
  Execute function 
-----------------------------------------------------------------*/
void FMR::execute()
{
  // ** Read input file to get run specs first ** //
  input->read_input_file(); 

  // ** Get atom information ** //
  if (input->read_restart) {
    // ** Load data from restart file ** // 
    input->read_restart_file();
  } else {
    // ** Load data from file ** //
    input->read_atoms_file();
  }

  // ** Run the requested calculation ** //
  run->run_calculation();

  // ** Final printing ** //
  if (master_rank) {
    printf("\n----------------------------------------------\n");
    printf("FMO-MS-RMD calculation completed successfully.\n");
    printf("----------------------------------------------\n");
  }
}

