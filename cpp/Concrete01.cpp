#include "Concrete01.h"
#include <string.h>
#include <math.h>
#include <float.h>

Concrete01::Concrete01
(int tag, double FPC, double EPSC0, double FPCU, double EPSCU):
   fpc(FPC), epsc0(EPSC0), fpcu(FPCU), epscu(EPSCU), 
   CminStrain(0.0), CendStrain(0.0),
   Cstrain(0.0), Cstress(0.0) 
{
  // Make all concrete parameters negative
  if (fpc > 0.0)
    fpc = -fpc;
  
  if (epsc0 > 0.0)
    epsc0 = -epsc0;
  
  if (fpcu > 0.0)
    fpcu = -fpcu;
  
  if (epscu > 0.0)
    epscu = -epscu;
  
  // Initial tangent
  double Ec0 = 2*fpc/epsc0;
  Ctangent = Ec0;
  CunloadSlope = Ec0;
  Ttangent = Ec0;
  
  // Set trial values
  this->revertToLastCommit();
}

Concrete01::~Concrete01(void)
{
  // Does nothing
}


int Concrete01::setTrialStrain (double strain, double strainRate)
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   TunloadSlope = CunloadSlope;
   Tstress = Cstress;
   Ttangent = Ctangent;
   Tstrain = Cstrain;

  // Determine change in strain from last converged state
  double dStrain = strain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON)
    return 0;

  // Set trial strain
  Tstrain = strain;
  
  // check for a quick return
  if (Tstrain > 0.0) {
    Tstress = 0;
    Ttangent = 0;
    return 0;
  }
  
  // Calculate the trial state given the change in strain
  // determineTrialState (dStrain);
  TunloadSlope = CunloadSlope;
  
  double tempStress = Cstress + TunloadSlope*Tstrain - TunloadSlope*Cstrain;
  
  // Material goes further into compression
  if (strain < Cstrain) {
    TminStrain = CminStrain;
    TendStrain = CendStrain;
    
    reload ();
    
    if (tempStress > Tstress) {
      Tstress = tempStress;
      Ttangent = TunloadSlope;
    }
  }
  
  // Material goes TOWARD tension
  else if (tempStress <= 0.0) {
    Tstress = tempStress;
    Ttangent = TunloadSlope;
  }
  
  // Made it into tension
  else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
  
  return 0;
}

void Concrete01::reload ()
{
  if (Tstrain <= TminStrain) {
    
    TminStrain = Tstrain;
    
    // Determine point on envelope
    envelope ();
    
    unload ();
  }
  else if (Tstrain <= TendStrain) {
    Ttangent = TunloadSlope;
    Tstress = Ttangent*(Tstrain-TendStrain);
  }
  else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
}

void Concrete01::envelope ()
{
  if (Tstrain > epsc0) {
    double eta = Tstrain/epsc0;
    Tstress = fpc*(2*eta-eta*eta);
    double Ec0 = 2.0*fpc/epsc0;
    Ttangent = Ec0*(1.0-eta);
  }
  else if (Tstrain > epscu) {
    Ttangent = (fpc-fpcu)/(epsc0-epscu);
    Tstress = fpc + Ttangent*(Tstrain-epsc0);
  }
  else {
    Tstress = fpcu;
    Ttangent = 0.0;
  }
}

void Concrete01::unload ()
{
  double tempStrain = TminStrain;
  
  if (tempStrain < epscu)
    tempStrain = epscu;
  
  double eta = tempStrain/epsc0;
  
  double ratio = 0.707*(eta-2.0) + 0.834;
  
  if (eta < 2.0)
    ratio = 0.145*eta*eta + 0.13*eta;
  
  TendStrain = ratio*epsc0;
  
  double temp1 = TminStrain - TendStrain;
  
  double Ec0 = 2.0*fpc/epsc0;
  
  double temp2 = Tstress/Ec0;
  
  if (temp1 > -DBL_EPSILON) {	// temp1 should always be negative
    TunloadSlope = Ec0;
  }
  else if (temp1 <= temp2) {
    TendStrain = TminStrain - temp1;
    TunloadSlope = Tstress/temp1;
  }
  else {
    TendStrain = TminStrain - temp2;
    TunloadSlope = Ec0;
  }
}

double Concrete01::getStress ()
{
   return Tstress;
}

double Concrete01::getStrain ()
{
   return Tstrain;
}

double Concrete01::getTangent ()
{
   return Ttangent;
}

int Concrete01::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CunloadSlope = TunloadSlope;
   CendStrain = TendStrain;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int Concrete01::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   TunloadSlope = CunloadSlope;

   // Recompute trial stress and tangent
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}