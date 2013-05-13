/* AWGL */
#ifndef FMR_POINTERS_H
#define FMR_POINTERS_H

#include "fmr.h"

namespace FMR_NS {

class Pointers {
 public:
  Pointers(FMR *ptr) :
   fmr(ptr), 
   input(ptr->input),
   atom(ptr->atom),
   run(ptr->run),
   matrix(ptr->matrix),
   state(ptr->state),
   dynamics(ptr->dynamics) {}
  virtual ~Pointers() {}

 protected:
  FMR *fmr;
  Atom *&atom;
  Input *&input;
  Run *&run;
  Matrix *&matrix;
  State *&state;
  Dynamics *&dynamics;

};

}

#endif
