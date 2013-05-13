/* AWGL */
#ifndef FMR_DYNAMICS_H
#define FMR_DYNAMICS_H

// Read options
#define DYNAMICS_ZERO  0 
#define DYNAMICS_READ  1
#define DYNAMICS_RAND  2 
#define DYNAMICS_MIN   0
#define DYNAMICS_MAX   2

// Thermostat options
#define THERMOSTAT_OFF        0
#define THERMOSTAT_LANGEVIN   1
#define THERMOSTAT_BERENDSEN  2

#include "pointers.h"

namespace FMR_NS {

class Dynamics : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Dynamics(FMR *);
   ~Dynamics();

   // ** Variables ** //
   int    nTimeSteps;             // Number of time steps to take 
   int    iCurrentStep;           // The current step we are on
   double dt;                     // Time step increment [atomic units]
   int    vOption;                // Option for starting velocity
   int    iThermostat;            // Thermostat option selected
   double targetTemperature;      // Target temperature of the run (Kelvin)
   double currentTemperature;     // The current instantaneous temperature (Kelvin)
   double EPotential;             // Total potential energy
   double EKinetic;               // Total nuclear kinetic energy (does not include electron kinetic energy)
   double ETotal;                 // EPotential + EKinetic
   int    vseed;                  // RNG seed for random velocities
   char   trajFile[256];          // File to write trajectory to
   double tau;                    // Thermostat time constant, used in Berendsen and Langevin

   // ** Functions ** //
   void   init();                   // Initializer function, calls others to set up MD run
   void   setVelocities();          // Set velocities based on vOption
   double computeKE();              // Compute the kinetic energy
   double computeTemperature();     // Compute the temperature
   void   runMolecularDynamics();   // Runs MD
   void   stepPosition();           // Move positions forward one step 
   void   halfStepVelocity(int);    // Propagate the velocities forward half a step
   void   stepReport(int);          // Print out stuff for time step
   void   writeTrajCoords(int);     // Write the current coordinates to trajectory file
   void   removeTransRotFromVelocity();  // Remove translation and rotation of atoms from velocity
   void   removeTransRotFromGradient();  // Remove translation and rotation of atoms from gradient
   void   applyThermostat();        // General function for applying thermostat for NVT dynamics
   void   thermostatLangevin();     // Langevin thermostat
   void   thermostatBerendsen();    // Berendsen thermostat

};

}

#endif
