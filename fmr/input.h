/* AWGL */
#ifndef FMR_INPUT_H
#define FMR_INPUT_H

#include "pointers.h"

namespace FMR_NS {

class Input : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Input(FMR *);
   ~Input();

   // ** Variables ** //
   char input_file[256];
   char atoms_file[256];
   char restart_file[256];
   int  read_restart;

   // ** Functions ** //
   void read_input_file();
   void read_atoms_file();
   void read_restart_file();
   void write_restart_file(); 

};

}

#endif
