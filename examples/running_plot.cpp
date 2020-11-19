#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "ckm.h"
#include "sm.h"
#include "pmns.h"
#include "rundown.cpp"

// convert g_x to alpha_x
double alpha(const double g){
  return g*g/4/M_PI;
    }

// write sm data to std output
struct write_out{
  void operator()(const sm &x, double t ){

    // calculate quark Yukawas
    ckm quarks(x.Yu, x.Yd);
    quarks.calculate();

    // calculate lepton Yukawas
    yukawa zero;
    pmns leptons(zero, x.Ye);
    leptons.calculate();
    
    // write to standard output
    std::cout << std::setprecision(10) << exp(t) << "   "
	      << alpha(x.g[0]) << "   "
      	      << alpha(x.g[1]) << "   "
      	      << alpha(x.g[2]) << "   "
	      << alpha(quarks.get_upyukawas()[2]) << "   "
	      << alpha(quarks.get_downyukawas()[2]) << "   "
	      << alpha(leptons.get_elyukawas()[2]) << std::endl;
  }
};

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
  
  // read the number of loops and set default to 2
  int nloops = (argc >= 2) ? atoi(argv[1]) : 1;
  std::cout << "#number of loops: " << nloops << std::endl;
  std::cout << "# alpha_1 alpha_2 alpha_3 alpha_t alpha_b alpha_tau " << std::endl;
  
  // define variables
  double MZ(91.1876);                      // EW scale
  double MGUT(1.e17);                      // GUT scale
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
  sm values(g,lambda,Yu,Yd,Ye,nloops);
  
  // set up stepper
  using namespace boost::numeric::odeint;
  typedef runge_kutta_fehlberg78< sm , double , sm , double , vector_space_algebra > stepper;

  // integrate RGEs from MZ to MGUT, note the observer function in the last argument of integrate_const
  int steps = integrate_const(stepper(), sm(), values , log(MZ) , log(MGUT), 0.1, write_out() );
          
  return 0;
}
