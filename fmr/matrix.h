/* AWGL */
#ifndef FMR_MATRIX_H
#define FMR_MATRIX_H

#include "pointers.h"
#include "math.h"

namespace FMR_NS {

// ** Inline switching functions ** //
inline double Sr(double r2, double rin2, double rout2) 
{
  // Switching function
  if (r2 < rin2) { 
    return 1.0;
  } else if (r2 > rout2) {
    return 0.0;
  } else {
    double denom = (rout2 - rin2);
    denom *= denom * denom;
    double tmp = rout2 - r2;
    return tmp*tmp * (rout2 + 2.0*r2 - 3.0*rin2) / denom;
  }
}

inline double dSrdr(double r2, double rin2, double rout2) 
{
  // Switching function derivative w/r r
  if (r2 < rin2) { 
    return 0.0;
  } else if (r2 > rout2) {
    return 0.0;
  } else {
    double denom = (rout2 - rin2);
    denom *= denom * denom;
    double tmp = rout2 - r2;
    double r = sqrt(r2);
    return -12.0*r * (rout2 - r2) * (r2 - rin2) / denom;
  }
}


// ** The class ** //
class Matrix : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Matrix(FMR *);
   ~Matrix();

   // ** Variables ** //
   double **H;          // The Hamiltonian matrix
   double **Evecs;      // Eigenvectors of H
   double *Evals;       // Eigenvalues of H
   double *Repulsion;   // Repulsion energy for each state
   double GSEnergy;     // Ground state energy
   double *GSCoeffs;    // Ground state eigenvector coefficients, points to Evecs[0]
   double ***HX;        // The gradient of the Hamiltonian matrix
   double **RepulsionX; // The gradient of the repulsion vector

   // ** Functions ** //
   void   buildH();                        // Builds the H matrix
   void   diagonalizeH();                  // Diagonalizes the current H matrix
   void   ComputeRepulsions();             // Fills up Repulsion array
   double ComputeRepulsionForState(int);   // Computes the repulsion energy for a given state
   void   ComputeCouplings();              // Fills up the off-diagonal couplings in H
   double ComputeCoupling(int,int);        // Compute coupling between I and J
   void   buildHX();                       // Builds the HX matrix
   void   ComputeRepulsionsX();            // Fills up RepulsionX array
   void   ComputeRepulsionXForState(int);  // Computes the repulsion gradient for a given state
   void   ComputeCouplingsX();             // Fills up the off-diagonal couplings in HX
   void   ComputeCouplingX(int,int);       // Compute coupling gradient between I and J
   void   ComputeHellmanFeynmanGradient(); // Computes ground state gradient via Hellman-Feynman
   double ComputeInternalPenaltyForState(int);   // Computes the internal penalty energy for a given state
   void   ComputeInternalPenaltyXForState(int);  // Computes the internal penalty gradient for a given state

};

}

#endif
