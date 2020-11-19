#include "pmns.h"
#include <cmath>

using namespace Eigen;
using namespace std;

void pmns::calculate(){  // core routine for calculating all parameters 
  
  // NOTE: eigen uses the same convention for rotation matrices as REAP/MPT
  
  // diagonalise charged leptons
  SelfAdjointEigenSolver<yukawa> EL(Ye.adjoint()*Ye);
  elyuk = EL.eigenvalues().array().sqrt().matrix();
  
  // diagonalise neutinos
  SelfAdjointEigenSolver<yukawa> NL(EL.eigenvectors().adjoint()*M.adjoint()*M*EL.eigenvectors());
  numasses = NL.eigenvalues().array().sqrt().matrix();
  
  // define PMNS matrix
  PMNS = NL.eigenvectors();
  
  // extract mixing parameters, MPT conventions
  double pi = 3.14159265358979323846;
  
  PMNSparameters[0] =  pi/2.; // default value for theta12
  if(PMNS(0,0).real() != 0.) PMNSparameters[0] = atan(abs(PMNS(0,1))/abs(PMNS(0,0)));
  
  PMNSparameters[1] = asin(abs(PMNS(0,2)));    // theta 13
  
  PMNSparameters[2] =  pi/2.; // default value for theta23
  if(PMNS(2,2).real() != 0.) PMNSparameters[2] = atan(abs(PMNS(1,2))/abs(PMNS(2,2)));
  
  PMNSparameters[3] = - arg(
			    (conj(PMNS(0,0))*PMNS(0,2)*PMNS(2,0)*conj(PMNS(2,2))
			     /cos(PMNSparameters[0])/cos(PMNSparameters[1])
			     /cos(PMNSparameters[1])/cos(PMNSparameters[2])
			     /sin(PMNSparameters[1])
			     + cos(PMNSparameters[0])*cos(PMNSparameters[2])
			     *sin(PMNSparameters[1]) )
			    / sin(PMNSparameters[0])/sin(PMNSparameters[2])
			    ); // Dirac phase should be positive in this convention!

  // always returns a positiv phase!
  if (PMNSparameters[3] < 0) PMNSparameters[3] = PMNSparameters[3]+2*pi;

  // calculate beta decay mass for Troitsk bound
  Vector3d PMNSabs2 = PMNS.row(0).cwiseAbs2();
  betadecaymass = sqrt(PMNSabs2.adjoint()*numasses.cwiseAbs2());
}

yukawa pmns::get_PMNS(){
  if (PMNS(0,0).real() == 0.) calculate();
  return PMNS;
}

Vector3d pmns::get_elyukawas(){
  if (elyuk[2]==0.) calculate();
  return elyuk;
}

Vector3d pmns::get_numasses(){
  if (numasses[2]==0.) calculate();
  return numasses;
}

Vector4d pmns::get_PMNSparameters(){
  if (PMNSparameters[0]==0.) calculate();
  return PMNSparameters;
}

Vector3d pmns::get_elmasses(){       // no argument: one Higgs doublet
  if (elyuk[0]==0.) calculate();
  return vev*elyuk;
}

Vector3d pmns::get_elmasses(double tanb){  // one argument: two Higgs doublets
  if (elyuk[0]==0.) calculate();
  return cos(atan(tanb))*vev*elyuk;
}

double pmns::get_betadecaymass(){
  if (numasses[2]==0.) calculate();
  return betadecaymass;
}

