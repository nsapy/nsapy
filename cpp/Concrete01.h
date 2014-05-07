// Description: This file contains the class definition for 
// Concrete01.h adapted from Concr1.f90 (Filippou)
//   - Modified Kent-Park envelope
//   - No tension
//   - Linear unloading/reloading
#include <float.h>

class Concrete01
{
 public:
  Concrete01 (int tag, double fpc, double eco, double fpcu, double ecu);
  ~Concrete01();

  const char *getClassType(void) const {return "Concrete01";};
  
  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain(void);      
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) {return 2.0*fpc/epsc0;}

  int commitState();
  int revertToLastCommit();

 protected:

 private:
  /*** Material Properties ***/
  double fpc;    // Compressive strength
  double epsc0;  // Strain at compressive strength
  double fpcu;   // Crushing strength
  double epscu;  // Strain at crushing strength
  
  /*** CONVERGED History Variables ***/
  double CminStrain;   // Smallest previous concrete strain (compression)
  double CunloadSlope; // Unloading (reloading) slope from CminStrain
  double CendStrain;   // Strain at the end of unloading from CminStrain
  
  /*** CONVERGED State Variables ***/
  double Cstrain;
  double Cstress;   
  double Ctangent;	// Don't need Ctangent other than for revert and sendSelf/recvSelf
  // Storing it is better than recomputing it!!!
  
  /*** TRIAL History Variables ***/
  double TminStrain;
  double TunloadSlope;
  double TendStrain;
  
  /*** TRIAL State Variables ***/
  double Tstrain;
  double Tstress;
  double Ttangent; // Not really a state variable, but declared here
  // for convenience
  
  void reload();
  void unload();
  void envelope();
};



