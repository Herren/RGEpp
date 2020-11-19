#include "nsm.h"
#include "ckm.h"
#include "pmns.h"

#define loopfactor 0.006332573977646111 // 1/(4 Pi)^2

using namespace Eigen;

// beta functions for sm
void nsm::operator()(const nsm &X, nsm &dX, const double){
  if (check()){

    // variables for powers of parameters for nloops >= 1
    gauge<3> g2(X.g.array().square().matrix());
    gauge<3> g3(X.g.array().cube().matrix());
    gauge<3> g4(g2.array().square().matrix());
    yukawa Yu2 = X.Yu.adjoint()*X.Yu; std::complex<double> Yu2Tr = Yu2.trace();
    yukawa Yd2 = X.Yd.adjoint()*X.Yd; std::complex<double> Yd2Tr = Yd2.trace();
    yukawa Ye2 = X.Ye.adjoint()*X.Ye; std::complex<double> Ye2Tr = Ye2.trace();
    yukawa Yu4 = Yu2*Yu2;             std::complex<double> Yu4Tr = Yu4.trace();
    yukawa Yd4 = Yd2*Yd2;             std::complex<double> Yd4Tr = Yd4.trace();
    yukawa Ye4 = Ye2*Ye2;             std::complex<double> Ye4Tr = Ye4.trace();
    yukawa Yu2Yd2 = Yu2*Yd2;          std::complex<double> Yu2Yd2Tr = Yu2Yd2.trace();
    yukawa Yd2Yu2 = Yd2*Yu2;

    dX.g[0] = loopfactor*((41./10.)*g3[0]);
    dX.g[1] = loopfactor*((-19./6.)*g3[1]);
    dX.g[2] = loopfactor*((-7.)*g3[2]);
    dX.Yu = loopfactor*((-3./2.)*X.Yu*Yd2 + (3./2.)*X.Yu*Yu2 + (3.)*X.Yu*Yd2Tr + X.Yu*Ye2Tr + (3.)*X.Yu*Yu2Tr + (-17./20.)*X.Yu*g2[0] + (-9./4.)*X.Yu*g2[1] + (-8.)*X.Yu*g2[2]);
    dX.Yd = loopfactor*((3./2.)*X.Yd*Yd2 + (-3./2.)*X.Yd*Yu2 + (3.)*X.Yd*Yd2Tr + X.Yd*Ye2Tr + (3.)*X.Yd*Yu2Tr + (-1./4.)*X.Yd*g2[0] + (-9./4.)*X.Yd*g2[1] + (-8.)*X.Yd*g2[2]);
    dX.Ye = loopfactor*((3./2.)*X.Ye*Ye2 + (3.)*X.Ye*Yd2Tr + X.Ye*Ye2Tr + (3.)*X.Ye*Yu2Tr + (-9./4.)*X.Ye*g2[0] + (-9./4.)*X.Ye*g2[1]);

    // NEW CODE for running of quartic couplings
    self<3> La2 = X.La.cwiseAbs2();
    
    dX.La[0] = loopfactor*((12.)*X.La[0]*Yd2Tr + (-24.)*Yd4Tr + (4.)*X.La[0]*Ye2Tr + (-8.)*Ye4Tr + (12.)*X.La[0]*Yu2Tr + (-24.)*Yu4Tr + (-9./5.)*X.La[0]*g2[0] + (-9.)*X.La[0]*g2[1] + (9./5.)*g2[0]*g2[1] + (27./50.)*g4[0] + (9./2.)*g4[1]
			   + 6.*La2[0] + 6.*La2[1]);
    dX.La[1] = loopfactor*( 4.*X.La[0]*X.La[1] + 4.*X.La[1]*X.La[2] );
    dX.La[2] = loopfactor*( 6.*La2[2] + 6.*La2[1]);


  } else {
    dX.setZero();
  }
};
