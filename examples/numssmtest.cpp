#include <iostream>
#include <cmath>
#include "numssm.h"
#include "ckm.h"
#include "pmns.h"
#include "rundown.cpp"

int main(int, char**){

  using namespace Eigen;

  // initiate parameters
  yukawa<3,3> Y10,Y126,zero;
  gauge<3> g;
  g << 0.718611, 0.71883, 0.719284; // defined at GUTscale

  double MGUT = 2e16;
  double MZ = 91.1876;
    
  // Dueck Rohdejohann best fit point for tan beta 50
  double tanb(50.);
  double r(1.87979);
  double rR(-1.09758e-15);
  std::complex<double> s(0.245295, 0.0935775);
  std::complex<double> H1(3.61945e-5, 0.);
  std::complex<double> H2(2.77898e-3, 0.);
  std::complex<double> H3(0.627274, 0.);
  std::complex<double> F1(4.13796e-6, 8.08833e-6);
  std::complex<double> F2(5.12296e-4, -5.07815e-4);
  std::complex<double> F3(-1.6274e-3, -2.47433e-3);
  std::complex<double> F4(-8.8214e-3, -5.30048e-5);
  std::complex<double> F5(-0.0170155, -0.019725);
  std::complex<double> F6(4.46e-3, -0.0583555);
  
  // construct matrices
  Y10 << H1,0,0, 0,H2,0, 0,0,H3;
  Y126 << F1,F2,F3, F2,F4,F5, F3,F5,F6;
  zero << 0,0,0,0,0,0,0,0,0;

  std::cout << Y10 << std::endl;
  std::cout << Y126 << std::endl;
      
  // initiate class for RGE evolution
  numssm values(g,                   // g
		r*(Y10+s*Y126),      // Yu
		Y10+Y126,            // Yd
		Y10-3.*Y126,         // Ye
		r*(Y10-3.*s*Y126),   // Yn
		zero,                // Ka
		Y126/rR);            // Mn
  
  // perform RGEs from M_Z to scale
  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78< numssm , double , numssm , double , vector_space_algebra > stepper;

  Vector3d seesaw = values.logthresholds();

  std::cout << seesaw.transpose() << std::endl;
  
  rundown<stepper,numssm,void> foo(values, log(MGUT), log(MZ));
  
  // extract quark masses and mixing
  class ckm quarks(values.Yu, values.Yd);
  quarks.calculate();
  
  // extract lepton masses and mixing
  yukawa<3,3> ML = values.getML(174.104, tanb);
			   
  class pmns leptons(ML, values.Ye);
  leptons.calculate();
  
  // write results to standard output
  std::cout << std::setprecision(5) << std::endl
	    << "SM parameters at MZ:" << std::endl
	    << "gauge couplings: " << values.g.transpose() << std::endl
    	    << "ckm parameters: " << quarks.get_CKMparameters().transpose() << std::endl
	    << "up-type masses: " << quarks.get_upmasses(tanb).transpose() << std::endl
	    << "down-type masses: " << quarks.get_downmasses(tanb).transpose() << std::endl
	    << "PMNS parameters: " << leptons.get_PMNSparameters().transpose() << std::endl
	    << "neutrino masses: " << leptons.get_numasses().transpose() << std:: endl
    	    << "charged lepton massess: " << leptons.get_elmasses(tanb).transpose() << std::endl
	    << std::endl;

  return(0);
}
