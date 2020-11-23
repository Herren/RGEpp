#include "eigentypes.h"
#include <Eigen/Eigenvalues>

#ifndef CKM_H
#define CKM_H

using namespace Eigen;

class ckm {
  
 private:
  yukawa Yu, Yd, CKM;             // Yukawa matrices
  Vector3d upyuk, downyuk;        // eigenvalues
  Vector4d CKMparameters;         // mixing parameters
  const double vev;               // convention "for fermion masses"
  
 public: 
  // constructor without arguments
 ckm() : 
    Yu(), Yd(),                // Yukawa matrices
    CKM(),                     // CKM matrix
    upyuk(), downyuk(),        // eigenvalues
    CKMparameters(),           // CKM parameters according to REAP/PT
    vev(174.104) {}            // convention "for fermion masses"
  
  // constructor with Yu and Yd as input parameter
 ckm(const yukawa Yuin, const yukawa Ydin) :
    Yu(Yuin), Yd(Ydin),                      // Yukawa matrices as input parameters
      CKM(), upyuk(),
      downyuk(),
      CKMparameters(),
      vev(174.104) {}

  // constructor with Yu, Yd and the vev as input parameter
 ckm(const yukawa Yuin, const yukawa Ydin, const double vevin) :
    Yu(Yuin), Yd(Ydin),                      // Yukawa matrices as input parameters
      CKM(), upyuk(),
      downyuk(),
      CKMparameters(),
      vev(vevin) {}
    
  void calculate();                // core routine

  // various member functions to return :
  // the CKM matrix
  // the eigenvalues of the Yukawas
  // the quark mixing parameters
  // quark masses (for one and two Higgs doublet models)
  
  yukawa get_CKM();
  Vector3d get_upyukawas();
  Vector3d get_downyukawas();
  Vector4d get_CKMparameters();
  Vector3d get_upmasses();
  Vector3d get_downmasses();
  Vector3d get_upmasses(const double vev);
  Vector3d get_downmasses(const double vev);
  
};

#endif
