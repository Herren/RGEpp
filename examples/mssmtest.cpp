#include <iostream>
#include <cmath>
#include "mssm.h"
#include "ckm.h"
#include "pmns.h"


// ckm matrix in the convention of the PDG
yukawa<3,3> ckm(double th12, double th13, double th23, double phi){
  yukawa<3,3> res;
  std::complex<double> iphi(0.,phi);
  res << cos(th12)*cos(th13), sin(th12)*cos(th13), sin(th13)*exp(-iphi), // first row
    -sin(th12)*cos(th23)-cos(th12)*sin(th23)*sin(th13)*exp(iphi),        // (2,1) element
    cos(th12)*cos(th23)-sin(th12)*sin(th23)*sin(th13)*exp(iphi),         // (2,2) element
    sin(th23)*cos(th13),                                                 // (2,3) element
    sin(th12)*sin(th23)-cos(th12)*cos(th23)*sin(th13)*exp(iphi),         // (3,1) element
    -cos(th12)*sin(th23)-sin(th12)*cos(th23)*sin(th13)*exp(iphi),        // (3,2) element
    cos(th23)*cos(th13);                                                 // (3,3) element
  return res;
}

int main(int, char**){

  using namespace Eigen;

  double MZ(91.1876);                      // value of MZ
  double tanb(10.);
  
  gauge<3> g_MZ;
  g_MZ <<  0.461425, 0.65184, 1.21272;

  yukawa<3,3> Yu_MZ, Yd_MZ, Ye_MZ,zero;

  // numbers cf. arXiv:1811.02895, Table 3
  Yu_MZ << 7.80222e-6, 0,0,0, 0.00364562, 0,0,0, 0.989661;
  Yd_MZ << 0.0000166293, 0,0,0, 0.000310436, 0,0,0, 0.0164568;
  Ye_MZ << 2.794745e-6, 0,0,0, 5.899863e-4, 0,0,0, 1.002950e-2;

  double th12(0.227035);                   // ckm angle theta 12
  double th13(0.00371224);                 // ckm angle theta 13
  double th23(0.0418112);                  // ckm angle theta 23
  double phi(1.14299);                     // ckm CP phase
    
  // rotate Yd
  Yd_MZ = Yd_MZ*ckm(th12,th13,th23,phi).adjoint();

  mssm values(g_MZ,Yu_MZ/sin(atan(tanb)),Yd_MZ/cos(atan(tanb)),Ye_MZ/cos(atan(tanb)));
  std::cout << values.getNloops() << std::endl;

  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78< mssm , double , mssm , double , vector_space_algebra > stepper;

  // perform RGEs from M_Z to 3 TeV
  int steps = integrate_adaptive( make_controlled<stepper>( 1E-12 , 1E-12 ) , mssm() , values , log(91.) , log(3000.) , 0.01 );

  // extract quark masses and mixing
  class ckm quarks(values.Yu, values.Yd);
  quarks.calculate();

  // extract lepton masses and mixing
  class pmns leptons(zero, values.Ye);
  leptons.calculate();
  
  // write results to standard output
  std::cout << std::setprecision(5) << std::endl
	    << "MSSM parameters at 3 TeV:" << std::endl
	    << "gauge couplings: " << values.g.transpose() << std::endl
	    << "up-type yukawas: " << quarks.get_upyukawas().transpose() << std::endl
	    << "down-type yukawas: " << quarks.get_downyukawas().transpose() << std::endl
    	    << "charged lepton yukawas: " << leptons.get_elyukawas().transpose() << std::endl
	    << "ckm parameters: " << quarks.get_CKMparameters().transpose() << std::endl
	    << std::endl;

    // 10 TeV
  steps = integrate_adaptive( make_controlled<stepper>( 1E-12 , 1E-12 ) , mssm() , values , log(3000.) , log(1.e4) , 0.01 );

  // extract quark masses and mixing
  class ckm quarks2(values.Yu, values.Yd);
  quarks2.calculate();

  // extract lepton masses and mixing
  class pmns leptons2(zero, values.Ye);
  leptons2.calculate();
  
  // write results to standard output
  std::cout << std::setprecision(5) << std::endl
	    << "MSSM parameters at 10 TeV:" << std::endl
	    << "gauge couplings: " << values.g.transpose() << std::endl
	    << "up-type yukawas: " << quarks2.get_upyukawas().transpose() << std::endl
	    << "down-type yukawas: " << quarks2.get_downyukawas().transpose() << std::endl
    	    << "charged lepton yukawas: " << leptons2.get_elyukawas().transpose() << std::endl
	    << "ckm parameters: " << quarks2.get_CKMparameters().transpose() << std::endl
	    << std::endl;

  // 10^16 GeV
  steps = integrate_adaptive( make_controlled<stepper>( 1E-12 , 1E-12 ) , mssm() , values , log(1.e4) , log(1.e16) , 0.01 );

  // extract quark masses and mixing
  class ckm quarks3(values.Yu, values.Yd);
  quarks3.calculate();

  // extract lepton masses and mixing
  class pmns leptons3(zero, values.Ye);
  leptons3.calculate();
  
  // write results to standard output
  std::cout << std::setprecision(5) << std::endl
	    << "MSSM parameters at 10^16 GeV:" << std::endl
	    << "gauge couplings: " << values.g.transpose() << std::endl
	    << "up-type yukawas: " << quarks3.get_upyukawas().transpose() << std::endl
	    << "down-type yukawas: " << quarks3.get_downyukawas().transpose() << std::endl
    	    << "charged lepton yukawas: " << leptons3.get_elyukawas().transpose() << std::endl
	    << "ckm parameters: " << quarks3.get_CKMparameters().transpose() << std::endl
	    << std::endl;  
   
  return(0);
}
