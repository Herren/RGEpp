#include <iostream>
#include <cmath>
#include "numssm.h"
#include "ckm.h"
#include "pmns.h"
#include "rundown.cpp"

int main(int, char**){

  using namespace Eigen;

  // define beginning and end of integration
  double MGUT = 2.e16;
  double MZ = 91.1876;

  // initiate and set model parameters
  yukawa<3,3> zero, Y1, Y2;
  gauge<3> g;

  Y1 << 1.e-4, 0., 0.,  0., 2.e-2, 0.,  0., 0., 0.4;
  Y2 << 0., -2.e-3, 3.e-2,  -2.e-3, 0., 0.5,  3.e-2, 0., 0.5;
  g << 0.7, 0.7, 0.7;

  // construct class for RGE evolution
  numssm values(g,                  // g
		Y1 + Y2,            // Yu
		Y1 - Y2,            // Yd
		Y1 - Y2,            // Ye
		Y1 + Y2,            // Yn
		zero,               // Ka
		Y1*1.e14);          // Mn
  
  // set the ODEint stepper for numerical integration
  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78< numssm , double , numssm , double , vector_space_algebra > stepper;

  // perform the integration from MGUT to MZ using the rundown template
  // it will consecutively integrate out right-handed neutrinos at their respective mass scale
  rundown<stepper,numssm,void> foo(values, log(MGUT), log(MZ));

  // define tan beta
  double tanb(60.);
    
  // extract quark masses and mixing at MZ
  class ckm quarks(values.Yu, values.Yd);
  quarks.calculate();
  
  // extract lepton masses and mixing at MZ
  yukawa<3,3> ML = values.getML(174.104, tanb);          // left-handed neutrino mass matrix
  class pmns leptons(ML, values.Ye);
  leptons.calculate();
  
  // write results to standard output
  std::cout << std::setprecision(5) << std::endl
	    << "SM parameters at MZ:" << std::endl
	    << "gauge couplings: " << values.g.transpose() << std::endl
    	    << "ckm parameters: " << quarks.get_CKMparameters().transpose() << std::endl
	    << "up-type masses (GeV): " << quarks.get_upmasses(tanb).transpose() << std::endl
	    << "down-type masses (GeV): " << quarks.get_downmasses(tanb).transpose() << std::endl
	    << "PMNS parameters: " << leptons.get_PMNSparameters().transpose() << std::endl
	    << "neutrino masses  (eV): " << leptons.get_numasses().transpose() << std:: endl
    	    << "charged lepton massess  (GeV): " << leptons.get_elmasses(tanb).transpose() << std::endl
	    << std::endl;

  return(0);
}
