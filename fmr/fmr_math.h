/* AWGL */
#ifndef FMR_MATH_H
#define FMR_MATH_H

#include "pointers.h"

namespace FMR_NS {

// This class holds various math constants, math functions, and
// molecular mechanics definitions and parameters

class Math : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Math(FMR *);
   ~Math();

   // ** Variables ** //
   int    seed;           // Random number generator seed
   bool   hasBeenSeeded;  // Has the RNG been seeded?
   // Math constants
   double pi;
   double twopi;
   // Conversion factors
   double b2a;            // Bohrs -> Angstroms
   double a2b;            // Angstroms -> Bohrs
   double au2kcal;        // hartrees -> kcal/mol
   double kcal2au;        // kcal/mol -> hartrees
   double kB_eV;          // Boltzmann constant in eV/Kelvin
   double kB_au;          // Boltzmann constant in hartree/Kelvin
   double amu2au;         // Atomic mass unit -> a.u. mass
   double fem2au;         // Femtosecond -> a.u. time 

   // Parameters
   double AA, BB, CC;     // FMO-MS-RMD diagonal and off-diagonal parameters 
   double rb;             // bond length threshold
   
   // ** Functions ** //
   void   rng_seed(int);               // Seed the random number generator
   double rng();                       // Return a random number in range 0 to 1
   double rngInRange(double, double);  // Return a random number in range min to max
   // Interface to LAPACK matrix inversion
   void   invertMatrix(double *, int); // Inverts an NxN matrix


};

}

#endif
