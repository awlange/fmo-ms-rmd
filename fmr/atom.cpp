/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/

Atom::Atom(FMR *fmr) : Pointers(fmr)
{
  // Basic initializations
  natoms       = 0;
  nfragments   = 0;
  nstates      = 0;
  prev_nstates = 0;
  ireactive    = 0;
  coord        = NULL;
  force        = NULL;
  veloc        = NULL;
  mass         = NULL;
  symbol       = NULL;
  fragment     = NULL;
  reactive     = NULL;
  available    = NULL;
  hop          = NULL;
  environment  = NULL;
  totalMass    = 1.0;

  // MM charges
  qO_SPCE = -0.8476;
  qH_SPCE = -(qO_SPCE) * 0.5;
  qO_hydronium = -0.5;
  qH_hydronium = (1.0 - qO_hydronium) / 3.0;
}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/

Atom::~Atom()
{
  if (coord     != NULL) delete [] coord; 
  if (force     != NULL) delete [] force; 
  if (veloc     != NULL) delete [] veloc; 
  if (mass      != NULL) delete [] mass; 
  if (symbol    != NULL) delete [] symbol; 
  if (fragment  != NULL) delete [] fragment; 
  if (reactive  != NULL) delete [] reactive; 
  if (available != NULL) delete [] available; 
  if (hop       != NULL) delete [] hop; 
}


/*-----------------------------------------------------------------
  Set all the atom masses to a.u. 
-----------------------------------------------------------------*/
void Atom::setAtomMasses()
{
  // Also compute the total mass
  totalMass = 0.0;
  for (int i=0; i<natoms; ++i) {
    mass[i] = getAtomMass(i) * fmr->math->amu2au;
    totalMass += mass[i];
  }
}
