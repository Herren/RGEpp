#include "ckm.h"
#include <cmath>

using namespace Eigen;
using namespace std;

void ckm::calculate() {  // core routine for calculating all parameters 
  
  // NOTE: eigen uses the same convention for rotation matrices as in REAP/MPT documentation
  
  // diagonalise up quarks
  SelfAdjointEigenSolver<yukawa<3,3> > UL(Yu.adjoint()*Yu);
  upyuk = UL.eigenvalues().array().sqrt().matrix();
  
  // diagonalise down quarks
  SelfAdjointEigenSolver<yukawa<3,3> > DL(UL.eigenvectors().adjoint()*Yd.adjoint()*Yd*UL.eigenvectors());
  downyuk = DL.eigenvalues().array().sqrt().matrix();
  
  // define CKM matrix
  CKM = DL.eigenvectors();
  
  // extract mixing parameters, MPT conventions
  double pi = 3.14159265358979323846;
  
  CKMparameters[0] =  pi/2.; // default value for theta12
  if(CKM(0,0).real() != 0.) CKMparameters[0] = atan(abs(CKM(0,1))/abs(CKM(0,0)));
  
  CKMparameters[1] = asin(abs(CKM(0,2)));    // theta 13
  
  CKMparameters[2] =  pi/2.; // default value for theta23
  if(CKM(2,2).real() != 0.) CKMparameters[2] = atan(abs(CKM(1,2))/abs(CKM(2,2)));
  
  CKMparameters[3] = - arg(( conj(CKM(0,0))*CKM(0,2)*CKM(2,0)*conj(CKM(2,2)) / cos(CKMparameters[0])/cos(CKMparameters[1])/cos(CKMparameters[1])/cos(CKMparameters[2])/sin(CKMparameters[1]) + cos(CKMparameters[0])*cos(CKMparameters[2])*sin(CKMparameters[1]) ) / sin(CKMparameters[0])/sin(CKMparameters[2])); // Dirac phase should be positive in this convention!
}

yukawa<3,3> ckm::get_CKM() {
  if (CKM(0,0).real() == 0.) calculate();
  return CKM;
}

Vector3d ckm::get_upyukawas() {
  if (upyuk[2]==0.) calculate();
  return upyuk;
}

Vector3d ckm::get_downyukawas() {
  if (downyuk[2]==0.) calculate();
  return downyuk;
}

Vector4d ckm::get_CKMparameters() {
  if (CKMparameters[0]==0.) calculate();
  return CKMparameters;
}

Vector3d ckm::get_upmasses() {       // no argument: one Higgs doublet
  if (upyuk[0]==0.) calculate();
  return vev*upyuk;
}

Vector3d ckm::get_downmasses() {
  if (downyuk[0]==0.) calculate();
  return vev*downyuk;
}

Vector3d ckm::get_upmasses(const double tanb) {  // one argument: two Higgs doublets
  if (upyuk[0]==0.) calculate();
  return sin(atan(tanb))*vev*upyuk;
}

Vector3d ckm::get_downmasses(const double tanb) {
  if (downyuk[0]==0.) calculate();
  return cos(atan(tanb))*vev*downyuk;
}
