#include <iostream>
#include <iomanip>
#include <string>
#include "ckm.h"
#include "thdm.h"
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


int main(int argc, char* argv[]){
  int nloops = (argc == 2) ? atoi(argv[1]) : 2;  
  // define variables
  double tanb;
  double MZ(91.1876);                      // EW scale
  double scale(3000.);                     // high energy scale
  yukawa<3,3> Yu,Yd,Ye,zero;                    // Yukawa matrices
  gauge<3> g;                                 // vector for gauge couplings
  std::complex<double> lambda(125.09/174.104/2.0,0.);   // quartic Higgs couplings

  // numbers cf. arXiv:1811.02895, Table 3
  Yu << 7.80222e-6, 0,0,0, 0.00364562, 0,0,0, 0.989661;
  Yd << 0.0000166293, 0,0,0, 0.000310436, 0,0,0, 0.0164568;
  Ye << 2.794745e-6, 0,0,0, 5.899863e-4, 0,0,0, 1.002950e-2;
  g <<  0.461425, 0.65184, 1.21272;
  
  double th12(0.227035);                   // ckm angle theta 12
  double th13(0.00371224);                 // ckm angle theta 13
  double th23(0.0418112);                  // ckm angle theta 23
  double phi(1.14299);                     // ckm CP phase

  // Set all quartic couplings equal
  self<5> La;
  La[0] = lambda;
  La[1] = lambda;
  La[2] = lambda;
  La[3] = lambda;
  La[4] = std::complex<double>(0.5,0.5);

  // Prepare Yukawa for 2HDM-Type I  
  tanb = 0.9;
  yukawa<3,3> Yui = Yu/sin(atan(tanb));
  yukawa<3,3> Ydi = Yd/sin(atan(tanb));
  yukawa<3,3> Yei = Ye/sin(atan(tanb));

  // Prepare Yukawa for 2HDM-Type II
  tanb = 50.0;
  yukawa<3,3> Yuii = Yu/sin(atan(tanb));
  yukawa<3,3> Ydii = Yd/cos(atan(tanb));
  yukawa<3,3> Yeii = Ye/cos(atan(tanb));

  // Prepare Yukawa for 2HDM-Type X  
  tanb = 100.0;
  yukawa<3,3> Yux = Yu/sin(atan(tanb));
  yukawa<3,3> Ydx = Yd/sin(atan(tanb));
  yukawa<3,3> Yex = Ye/cos(atan(tanb));

  // Prepare Yukawa for 2HDM-Type Y
  tanb = 25.0;
  yukawa<3,3> Yuy = Yu/sin(atan(tanb));
  yukawa<3,3> Ydy = Yd/cos(atan(tanb));
  yukawa<3,3> Yey = Ye/sin(atan(tanb));

  // rotate Yd
  Ydi = Ydi*ckm(th12,th13,th23,phi).adjoint();
  Ydii = Ydii*ckm(th12,th13,th23,phi).adjoint();
  Ydx = Ydx*ckm(th12,th13,th23,phi).adjoint();
  Ydy = Ydy*ckm(th12,th13,th23,phi).adjoint();

  // add data to the RGE class
  thdmi valuesi(g,La,Yui,Ydi,Yei,nloops);
  thdmii valuesii(g,La,Yuii,Ydii,Yeii,nloops);
  thdmx valuesx(g,La,Yux,Ydx,Yex,nloops);
  thdmy valuesy(g,La,Yuy,Ydy,Yey,nloops);
  
  // perform RGEs from M_Z to scale
  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78< thdmi , double , thdmi , double , vector_space_algebra > stepperi;
  int stepsi = integrate_adaptive( make_controlled<stepperi>( 1E-12 , 1E-12 ) , thdmi() , valuesi , log(MZ) , log(scale) , 0.01 );
  typedef runge_kutta_fehlberg78< thdmii , double , thdmii , double , vector_space_algebra > stepperii;
  int stepsii = integrate_adaptive( make_controlled<stepperii>( 1E-12 , 1E-12 ) , thdmii() , valuesii , log(MZ) , log(scale) , 0.01 );
  typedef runge_kutta_fehlberg78< thdmx , double , thdmx , double , vector_space_algebra > stepperx;
  int stepsx = integrate_adaptive( make_controlled<stepperx>( 1E-12 , 1E-12 ) , thdmx() , valuesx , log(MZ) , log(scale) , 0.01 );
  typedef runge_kutta_fehlberg78< thdmy , double , thdmy , double , vector_space_algebra > steppery;
  int stepsy = integrate_adaptive( make_controlled<steppery>( 1E-12 , 1E-12 ) , thdmy() , valuesy , log(MZ) , log(scale) , 0.01 );

  // extract quark masses and mixing
  class ckm quarksi(valuesi.Yu, valuesi.Yd);
  class ckm quarksii(valuesii.Yu, valuesii.Yd);
  class ckm quarksx(valuesx.Yu, valuesx.Yd);
  class ckm quarksy(valuesy.Yu, valuesy.Yd);
  quarksi.calculate();
  quarksii.calculate();
  quarksx.calculate();
  quarksy.calculate();

  // extract lepton masses and mixing
  class pmns leptonsi(zero, valuesi.Ye);
  class pmns leptonsii(zero, valuesii.Ye);
  class pmns leptonsx(zero, valuesx.Ye);
  class pmns leptonsy(zero, valuesy.Ye);
  leptonsi.calculate();
  leptonsii.calculate();
  leptonsx.calculate();
  leptonsy.calculate();
  
  // write results to standard output
  std::cout << std::setprecision(5) << std::endl
	    << "2HDM Type I parameters at " << scale << " GeV:" << std::endl
	    << "gauge couplings: " << valuesi.g.transpose() << std::endl
	    << "up-type yukawas: " << quarksi.get_upyukawas().transpose() << std::endl
	    << "down-type yukawas: " << quarksi.get_downyukawas().transpose() << std::endl
    	    << "charged lepton yukawas: " << leptonsi.get_elyukawas().transpose() << std::endl
	    << "ckm parameters: " << quarksi.get_CKMparameters().transpose() << std::endl
	    << " Higgs selfcoupling: " << valuesi.La.transpose() << std::endl
	    << std::endl;

  std::cout << std::setprecision(5) << std::endl
	    << "2HDM Type II parameters at " << scale << " GeV:" << std::endl
	    << "gauge couplings: " << valuesii.g.transpose() << std::endl
	    << "up-type yukawas: " << quarksii.get_upyukawas().transpose() << std::endl
	    << "down-type yukawas: " << quarksii.get_downyukawas().transpose() << std::endl
    	    << "charged lepton yukawas: " << leptonsii.get_elyukawas().transpose() << std::endl
	    << "ckm parameters: " << quarksii.get_CKMparameters().transpose() << std::endl
	    << " Higgs selfcoupling: " << valuesii.La.transpose() << std::endl
	    << std::endl;

  std::cout << std::setprecision(5) << std::endl
	    << "2HDM Type X parameters at " << scale << " GeV:" << std::endl
	    << "gauge couplings: " << valuesx.g.transpose() << std::endl
	    << "up-type yukawas: " << quarksx.get_upyukawas().transpose() << std::endl
	    << "down-type yukawas: " << quarksx.get_downyukawas().transpose() << std::endl
    	    << "charged lepton yukawas: " << leptonsx.get_elyukawas().transpose() << std::endl
	    << "ckm parameters: " << quarksx.get_CKMparameters().transpose() << std::endl
	    << " Higgs selfcoupling: " << valuesx.La.transpose() << std::endl
	    << std::endl;

  std::cout << std::setprecision(5) << std::endl
	    << "2HDM Type Y parameters at " << scale << " GeV:" << std::endl
	    << "gauge couplings: " << valuesy.g.transpose() << std::endl
	    << "up-type yukawas: " << quarksy.get_upyukawas().transpose() << std::endl
	    << "down-type yukawas: " << quarksy.get_downyukawas().transpose() << std::endl
    	    << "charged lepton yukawas: " << leptonsy.get_elyukawas().transpose() << std::endl
	    << "ckm parameters: " << quarksy.get_CKMparameters().transpose() << std::endl
	    << " Higgs selfcoupling: " << valuesy.La.transpose() << std::endl
	    << std::endl;
      
  return 0;
}
