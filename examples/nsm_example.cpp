#include <iostream>
#include <iomanip>
#include <string>
#include "nsm.h"
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

  // define variables
  double MZ(91.1876);                      // EW scale
  yukawa<3,3> Yu,Yd,Ye,zero;                    // Yukawa matrices
  gauge<3> g;                              // vector for gauge couplings

  // numbers cf. arXiv:1811.02895, Table 3
  Yu << 7.80222e-6, 0,0,0, 0.00364562, 0,0,0, 0.989661;
  Yd << 0.0000166293, 0,0,0, 0.000310436, 0,0,0, 0.0164568;
  Ye << 2.794745e-6, 0,0,0, 5.899863e-4, 0,0,0, 1.002950e-2;
  g <<  0.461425, 0.65184, 1.21272;
  
  double th12(0.227035);                   // ckm angle theta 12
  double th13(0.00371224);                 // ckm angle theta 13
  double th23(0.0418112);                  // ckm angle theta 23
  double phi(1.14299);                     // ckm CP phase
  
  // rotate Yd
  Yd = Yd*ckm(th12,th13,th23,phi).adjoint();

  // NEW CODE
  std::complex<double> l1(125.09/174.104,0.);   // quartic doublet coupling
  std::complex<double> l2(1,0.);                // singlet-doublet coupling
  std::complex<double> l3(0.1,0.);              // quartic single coupling

  // add data to the RGE class
  nsm values(g,l1,l2,l3,Yu,Yd,Ye);
  
  // perform RGEs from M_Z to scale
  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78< nsm , double , nsm , double , vector_space_algebra > stepper;
  int steps = integrate_adaptive( make_controlled<stepper>( 1E-10 , 1E-10 ) , nsm() , values , log(MZ) , log(1000.) , 0.01 );

  
  // write results to standard output
  std::cout << std::setprecision(8) << std::endl
	    << "nSM parameters at " << " 1000 GeV:" << std::endl
	    << "Higgs selfcouplings: " << values.La.transpose() << std::endl
	    << std::endl;

      
  return 0;
}
