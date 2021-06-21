#include <iostream>
#include <iomanip>
#include <string>
#include "sm.h"
#include "ckm.h"
#include "pmns.h"

// ckm matrix in the convention of the PDG
yukawa ckm(double th12, double th13, double th23, double phi){
  yukawa res;
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


int main(int argc, char* argv[]){

  // read the final renormalisation scale from comand line and set the default to 3 TeV
  double scale = (argc >= 2) ? atof(argv[1]) : 3000.;

  // read the number of loops for the beta function and set default to 2
  int nloops = (argc >= 3) ? atoi(argv[2]) : 2;
  
  // read flag for L-(L-1)-(L-2) ordering and set default to false
  bool weylordering = (argc >= 4) ? atoi(argv[3]) : false;

  // define variables
  double MZ(91.1876);                      // EW scale
  yukawa Yu,Yd,Ye,zero;                    // Yukawa matrices
  gauge<3> g;                                 // vector for gauge couplings
  std::complex<double> lambda(125.09/174.104,0.);   // quartic Higgs coupling

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

  // add data to the RGE class
  sm values(g,lambda,Yu,Yd,Ye,nloops,weylordering);
  
  // perform RGEs from M_Z to scale
  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78< sm , double , sm , double , vector_space_algebra > stepper;
  int steps = integrate_adaptive( make_controlled<stepper>( 1E-10 , 1E-10 ) , sm() , values , log(MZ) , log(scale) , 0.01 );

  // extract quark masses and mixing
  class ckm quarks(values.Yu, values.Yd);
  quarks.calculate();

  // extract lepton masses and mixing
  class pmns leptons(zero, values.Ye);
  leptons.calculate();
  
  // write results to standard output
  std::cout << std::setprecision(8) << std::endl
	    << "SM parameters at " << scale << " GeV and " << values.getNloops() << " loops:" << std::endl
	    << "gauge couplings: " << values.g.transpose() << std::endl
	    << "up-type yukawas: " << quarks.get_upyukawas().transpose() << std::endl
	    << "down-type yukawas: " << quarks.get_downyukawas().transpose() << std::endl
    	    << "charged lepton yukawas: " << leptons.get_elyukawas().transpose() << std::endl
	    << "ckm parameters: " << quarks.get_CKMparameters().transpose() << std::endl
	    << "Higgs quartic coupling: " << values.La[0] << std::endl
	    << std::endl;

      
  return 0;
}
