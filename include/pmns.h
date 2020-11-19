#include "eigentypes.h"
#include <Eigen/Eigenvalues>

#ifndef PMNS_H
#define PMNS_H

using namespace Eigen;

class pmns{
  
 private:
  yukawa Ye, M, PMNS;             // Yukawa matrices
  Vector3d elyuk, numasses;       // eigenvalues
  Vector4d PMNSparameters;        // mixing parameters
  const double vev;               // convention "for fermion masses"
  double betadecaymass;           // for Troitsk bound
  
 public:
  
  // constructor without arguments
 pmns() : Ye(), M(),           // Yukawa matrices
    PMNS(),                    // CKM matrix
    elyuk(), numasses(),        // eigenvalues
    PMNSparameters(),          // CKM parameters according to REAP/PT
    vev(174.104) {}            // convention "for fermion masses"
  
  // constructor with Ye and M as input parameter
 pmns(yukawa Min, yukawa Yein) :
  Ye(Yein), M(Min),                      // Yukawa matrices as input parameters
    PMNS(), elyuk(), numasses(), PMNSparameters(), vev(174.104) {}

  // constructor with Ye, M and vev as input parameter
 pmns(yukawa Min, yukawa Yein, double vevin) :
  Ye(Yein), M(Min),                      // Yukawa matrices as input parameters
    PMNS(), elyuk(), numasses(), PMNSparameters(), vev(vevin) {}
  
  void calculate();                // core routine

  // various member functions to return :
  // the CKM matrix
  // the eigenvalues of the Yukawas
  // the quark mixing parameters
  // quark masses (for one and two Higgs doublet models)
  
  yukawa get_PMNS();
  Vector3d get_elyukawas();
  Vector3d get_numasses();
  Vector4d get_PMNSparameters();
  Vector3d get_elmasses();
  Vector3d get_elmasses(double vev);
  double get_betadecaymass();
  
};

#endif
