/* AWGL */
#ifndef FMR_ATOM_H
#define FMR_ATOM_H

#include "pointers.h"

#define MAX_STATES 30

namespace FMR_NS {

class Atom : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Atom(FMR *);
   ~Atom();

   // ** Variables ** //
   int natoms;           // Number of atoms (static)
   int nfragments;       // Number of fragments (static)
   int nstates;          // Number of fragmentation states (dynamic)
   int prev_nstates;     // Number of fragmentation states in previous step 
   int ireactive;        // Index of initial reactive fragment
   double *coord;        // Cartesian coordinates of atoms in Angstroms
   double *force;        // Forces on atoms (IMPORTANT: actually the gradient)
   double *veloc;        // Velocities of atoms
   double *mass;         // Mass of atoms
   double totalMass;     // Sum of all atom masses
   char *symbol;         // Atomic symbol of atoms
   int *fragment;        // Fragment index of atom for each fragmentation state
   int *reactive;        // Index specifying if atom belongs to a reactive fragment in each state
   int *available;       // Availability of atom in state search
   int *hop;             // How many hops this atom is from the pivot state reactive fragment
   int *environment;     // Do I belong to an environment fragment? 1=yes, 0=no
   // MM charge parameters
   double qO_SPCE;        // SPC/E charge oxygen 
   double qH_SPCE;        // SPC/E charge hydrogen
   double qO_hydronium;   // Hydronium charge oxygen
   double qH_hydronium;   // Hydronium charge hydrogen


   // ** Functions ** //
   bool AtomInFragment(int iatom, int ifrag, int istate){ 
     return fragment[istate*natoms + iatom] == ifrag;
   }

   double getCharge(int iatom, int istate) {
     double mmq = 0.0;
     if (symbol[iatom] == 'H') {
       if (reactive[istate*natoms + iatom]) mmq = qH_hydronium; 
       else                                 mmq = qH_SPCE;
     } else if (symbol[iatom] == 'O') {
       if (reactive[istate*natoms + iatom]) mmq = qO_hydronium; 
       else                                 mmq = qO_SPCE;
     } else {
       printf("Charge not found for atom %d of state %d\n", iatom, istate);
     }

     return mmq;
   }

   double getAtomMass(int iatom) {
     // Return atom mass in AMU. Must convert later to a.u. 
     double mass = 0.0;
     if      (symbol[iatom] == 'H') mass = 1.00783;
     else if (symbol[iatom] == 'O') mass = 15.99491;
     return mass;     
   }

   void setAtomMasses(); // in atom.cpp

};

}

#endif
