#ifndef Steel01_h
#define Steel01_h

// Default values for isotropic hardening parameters a1, a2, a3, and a4
#define STEEL_01_DEFAULT_A1        0.0
#define STEEL_01_DEFAULT_A2       55.0
#define STEEL_01_DEFAULT_A3        0.0
#define STEEL_01_DEFAULT_A4       55.0

class Steel01
{
  public:
    Steel01(int tag, double fy, double E0, double b,
       double a1 = STEEL_01_DEFAULT_A1, double a2 = STEEL_01_DEFAULT_A2,
       double a3 = STEEL_01_DEFAULT_A3, double a4 = STEEL_01_DEFAULT_A4);
    ~Steel01();

    const char *getClassType(void) const {return "Steel01";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return E0;};

    int commitState(void);
    int revertToLastCommit(void);    
    
 protected:
    
 private:
    /*** Material Properties ***/
    double fy;  // Yield stress
    double E0;  // Initial stiffness
    double b;   // Hardening ratio (b = Esh/E0)
    double a1;
    double a2;
    double a3;
    double a4;  // a1 through a4 are coefficients for isotropic hardening
    
    /*** CONVERGED History Variables ***/
    double CminStrain;  // Minimum strain in compression
    double CmaxStrain;  // Maximum strain in tension
    double CshiftP;     // Shift in hysteresis loop for positive loading
    double CshiftN;     // Shift in hysteresis loop for negative loading
    int Cloading;       // Flag for loading/unloading
                        // 1 = loading (positive strain increment)
                        // -1 = unloading (negative strain increment)
                        // 0 initially

    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cstress;
    double Ctangent;    

    /*** TRIAL History Variables ***/
    double TminStrain;
    double TmaxStrain;
    double TshiftP;
    double TshiftN;
    int Tloading;
    
    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience

    // Calculates the trial state variables based on the trial strain
    void determineTrialState (double dStrain);
};

#endif
