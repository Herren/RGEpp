#include <iostream>
#include <cmath>
#include "nusm.h"
#include "ckm.h"
#include "pmns.h"
#include "rundown.cpp"

struct yn_observer{
  void operator()(const nusm &x, double t){
    std::cout << t << "   " << x.Yn(1,1) << std::endl;
  }
};


int main(int, char**){

  using namespace Eigen;

  // initiate parameters
  std::complex<double> lambda(0.,0.);
  yukawa Y10,Y126,zero;

  gauge<3> g;
  g <<  0.579086, 0.521435, 0.526626; // gauge couplings at the GUT scale
  
  double MGUT = 2e16;
  double MZ = 91.1876;
  
  // Dueck Rohdejohann best fit point for MN, RGE
  // see arxiv:1306.4468
  double r(-63.9043);
  double rR(3.39236e-16);
  std::complex<double> s(0.409807, -0.0420522);
  std::complex<double> H1(1.15249e-6, 0.);
  std::complex<double> H2(6.71983e-5, 0.);
  std::complex<double> H3(6.70159e-3, 0.);
  std::complex<double> F1(-2.25817e-6, 7.40559e-6);
  std::complex<double> F2(1.22057e-5, -1.39062e-5);
  std::complex<double> F3(-1.49653e-4, +8.30809e-5);
  std::complex<double> F4(-2.06635e-4, -1.34686e-5);
  std::complex<double> F5(+3.76355e-4, +2.15049e-4);
  std::complex<double> F6(-7.01333e-4, -7.53673e-4);

  // construct matrices
  Y10 << H1,0,0, 0,H2,0, 0,0,H3;
  Y126 << F1,F2,F3, F2,F4,F5, F3,F5,F6;

  // initiate class for RGE evolution
  nusm values(g,                   // g
	      lambda,              // lambda
	      r*(Y10+s*Y126),      // Yu
	      Y10+Y126,            // Yd
	      Y10-3.*Y126,         // Ye
	      r*(Y10-3.*s*Y126),   // Yn
	      zero,                // Ka
	      Y126/rR);            // Mn

  

  Eigen::Vector3d seesaw = values.logthresholds();

  // perform RGEs from MGUT to MZ
  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78< nusm , double , nusm , double , vector_space_algebra > stepper;
  //  rundown<stepper, nusm, yn_observer> foo(values, log(MGUT), log());

  int steps = integrate_adaptive( make_controlled<stepper>( 1E-12 , 1E-12 ) , nusm() , values , log(MGUT) , seesaw[2] , -0.01 );

  values.integrate_out(2);
  
  steps = integrate_adaptive( make_controlled<stepper>( 1E-12 , 1E-12 ) , nusm() , values , seesaw[2] , seesaw[1] , -0.01 );
  
  values.integrate_out(1);
  
  steps = integrate_adaptive( make_controlled<stepper>( 1E-12 , 1E-12 ) , nusm() , values , seesaw[1] , seesaw[0] , -0.01 );
  
  values.integrate_out(0);
  
  steps = integrate_adaptive( make_controlled<stepper>( 1E-12 , 1E-12 ) , nusm() , values , seesaw[0] , log(MZ) , -0.01 );
  
  std::cout << "Kappa: " << values.Ka << std::endl;
  std::cout << "Mn: " << values.Mn << std::endl;
  std::cout << "Yn: " << values.Yn << std::endl;
  
  // extract quark masses and mixing
  class ckm quarks(values.Yu, values.Yd);
  quarks.calculate();
  
  // extract lepton masses and mixing
  yukawa ML =   -0.5e9*174.104*174.104*values.Ka;
			   
  class pmns leptons(ML, values.Ye);
  leptons.calculate();
  
  // write results to standard output
  std::cout << std::setprecision(5) << std::endl
	    << "SM parameters at MZ:" << std::endl
	    << "gauge couplings: " << values.g.transpose() << std::endl
    	    << "ckm parameters: " << quarks.get_CKMparameters().transpose() << std::endl
	    << "up-type masses: " << quarks.get_upmasses().transpose() << std::endl
	    << "down-type masses: " << quarks.get_downmasses().transpose() << std::endl
	    << "PMNS parameters: " << leptons.get_PMNSparameters().transpose() << std::endl
	    << "neutrino masses: " << leptons.get_numasses().transpose() << std:: endl
    	    << "charged lepton massess: " << leptons.get_elmasses().transpose() << std::endl
	    << std::endl;
   
  return(0);
}
