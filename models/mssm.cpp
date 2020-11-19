#include "mssm.h"

#define loopfactor 0.006332573977646111 // 1/(4 Pi)^2
#define loopfactor2 0.00004010149318236069 // 1/(4 Pi)^4
#define loopfactor3 0.00000025394567219137 // 1/(4 Pi)^6
#define z3 1.2020569031595942854 // Zeta(3)

using namespace Eigen;

// beta functions for mssm
void mssm::operator()(const mssm &X, mssm &dX, const double){
  if (check()){

    // variables for powers of parameters for nloops >= 1
    gauge<3> g2(X.g.array().square().matrix());
    gauge<3> g3(X.g.array().cube().matrix());
    yukawa Yu2 = X.Yu.adjoint()*X.Yu; yukawa Ytu2 = X.Yu*X.Yu.adjoint(); std::complex<double> Yu2Tr = Yu2.trace();
    yukawa Yd2 = X.Yd.adjoint()*X.Yd; yukawa Ytd2 = X.Yd*X.Yd.adjoint(); std::complex<double> Yd2Tr = Yd2.trace();
    yukawa Ye2 = X.Ye.adjoint()*X.Ye; yukawa Yte2 = X.Ye*X.Ye.adjoint(); std::complex<double> Ye2Tr = Ye2.trace();

    dX.g[0] = loopfactor*((33./5.)*g3[0]);
    dX.g[1] = loopfactor*(g3[1]);
    dX.g[2] = loopfactor*((-3.)*g3[2]);
    dX.Yu = loopfactor*((-13./15.)*g2[0]*X.Yu + (-3.)*g2[1]*X.Yu + (-16./3.)*g2[2]*X.Yu + (3.)*Yu2Tr*X.Yu + (3.)*X.Yu*Yu2 + (1.)*X.Yu*Yd2);
    dX.Yd = loopfactor*((-7./15.)*g2[0]*X.Yd + (-3.)*g2[1]*X.Yd + (-16./3.)*g2[2]*X.Yd + (3.)*Yd2Tr*X.Yd + (1.)*Ye2Tr*X.Yd + (1.)*X.Yd*Yu2 + (3.)*X.Yd*Yd2);
    dX.Ye = loopfactor*((-9./5.)*g2[0]*X.Ye + (-3.)*g2[1]*X.Ye + (1.)*Ye2Tr*X.Ye + (3.)*Yd2Tr*X.Ye + (3.)*X.Ye*Ye2);

    if(nloops > 1) {
      // variables for powers of parameters for nloops >= 2
      gauge<3> g4(g2.array().square().matrix());
      yukawa Yu4 = Yu2*Yu2; yukawa Ytu4 = Ytu2*Ytu2; std::complex<double> Yu4Tr = Yu4.trace();
      yukawa Yd4 = Yd2*Yd2; yukawa Ytd4 = Ytd2*Ytd2; std::complex<double> Yd4Tr = Yd4.trace();
      yukawa Ye4 = Ye2*Ye2; yukawa Yte4 = Yte2*Yte2; std::complex<double> Ye4Tr = Ye4.trace();
      yukawa Yu2Yd2 = Yu2*Yd2; std::complex<double> Yu2Yd2Tr = Yu2Yd2.trace();

      dX.g[0] += loopfactor2*g3[0]*((199./25.)*g2[0] + (27./5.)*g2[1] + (88./5.)*g2[2] + (-26./5.)*Yu2Tr + (-14./5.)*Yd2Tr + (-18./5.)*Ye2Tr).real();
      dX.g[1] += loopfactor2*g3[1]*((9./5.)*g2[0] + (25.)*g2[1] + (24.)*g2[2] + (-6.)*Yu2Tr + (-6.)*Yd2Tr + (-2.)*Ye2Tr).real();
      dX.g[2] += loopfactor2*g3[2]*((11./5.)*g2[0] + (9.)*g2[1] + (14.)*g2[2] + (-4.)*Yu2Tr + (-4.)*Yd2Tr).real();


      // small discrepancy w.r.t. old code in Yd?
      dX.Yu += loopfactor2*(((-2.)*Ytu4 + (-6.)*Ytu2*Yu2Tr + (-2.)*X.Yu*Yd2*X.Yu.adjoint() + (6.)*g2[1]*Ytu2 + (-2./5.)*g2[0]*Ytu2 + (856./225.)*g4[0] + (128./45.)*g2[0]*g2[2] + (-8./9.)*g4[2])*X.Yu
                         + X.Yu*((-2.)*Yu4 + (-3.)*Yu2*Yu2Tr + (-2.)*Yd4 + (-3.)*Yd2*Yd2Tr + (-1.)*Yd2*Ye2Tr + (4./5.)*g2[0]*Yu2 + (2./5.)*g2[0]*Yd2 + (1./10.)*g2[0]*g2[1] + (8./45.)*g2[0]*g2[2]
                                 + (8.)*g2[1]*g2[2] + (199./900.)*g4[0] + (15./4.)*g4[1] + (-8./9.)*g4[2])
                         + X.Yu*((-9.)*Yu4Tr + (-3.)*Yu2Yd2Tr + (16.)*g2[2]*Yu2Tr + (4./5.)*g2[0]*Yu2Tr + (9./10.)*g2[0]*g2[1] + (207./100.)*g4[0] + (15./4.)*g4[1]));
      dX.Yd += loopfactor2*(((-2.)*Ytd4 + (-6.)*Ytd2*Yd2Tr + (-2.)*X.Yd*Yu2*X.Yd.adjoint() + (-2.)*Ytd2*Ye2Tr + (6.)*g2[1]*Ytd2 + (2./5.)*g2[0]*Ytd2 + (202./225.)*g4[0] + (32./45.)*g2[0]*g2[2] + (-8./9.)*g4[2])*X.Yd
                         + X.Yd*((-2.)*Yu4 + (-3.)*Yu2*Yu2Tr + (-2.)*Yd4 + (-3.)*Yd2*Yd2Tr + (-1.)*Yd2*Ye2Tr + (4./5.)*g2[0]*Yu2 + (2./5.)*g2[0]*Yd2 + (1./10.)*g2[0]*g2[1] + (8./45.)*g2[0]*g2[2]
                                 + (8.)*g2[1]*g2[2] + (199./900.)*g4[0] + (15./4.)*g4[1] + (-8./9.)*g4[2])
                         + X.Yd*((-9.)*Yd4Tr + (-3.)*Yu2Yd2Tr + (-3.)*Ye4Tr + (16.)*g2[2]*Yd2Tr + (-2./5.)*g2[0]*Yd2Tr + (6./5.)*g2[0]*Ye2Tr + (9./10.)*g2[0]*g2[1] + (207./100.)*g4[0] + (15./4.)*g4[1]));
      dX.Ye += loopfactor2*(((-2.)*Yte4 + (-6.)*Yte2*Yd2Tr + (-2.)*Yte2*Ye2Tr + (6.)*g2[1]*Yte2 + (-6./5.)*g2[0]*Yte2 + (234./25.)*g4[0])*X.Ye
                         + X.Ye*((-2.)*Ye4 + (-3.)*Ye2*Yd2Tr + (-1.)*Ye2*Ye2Tr + (6./5.)*g2[0]*Ye2 + (9./10.)*g2[0]*g2[1] + (207./100.)*g4[0] + (15./4.)*g4[1])
                         + X.Ye*((-9.)*Yd4Tr + (-3.)*Yu2Yd2Tr + (-3.)*Ye4Tr + (16.)*g2[2]*Yd2Tr + (-2./5.)*g2[0]*Yd2Tr + (6./5.)*g2[0]*Ye2Tr + (9./10.)*g2[0]*g2[1] + (207./100.)*g4[0] + (15./4.)*g4[1]));
      if(nloops > 2) {
        // variables for powers of parameters for nloops >= 3
        dX.g[0] += loopfactor3*g3[0]*((84./5.)*Yu4Tr + (18.)*Yu2Tr*Yu2Tr + (58./5.)*Yu2Yd2Tr + (54./5.)*Yd4Tr + (36./5.)*Yd2Tr*Yd2Tr + (84./5.)*Ye2Tr*Yd2Tr + (54./5.)*Ye4Tr + (24./5.)*Ye2Tr*Ye2Tr
                                    + ((-352./15.)*g2[2] + (-87./5.)*g2[1] + (-169./75.)*g2[0])*Yu2Tr + ((-256./15.)*g2[2] + (-33./5.)*g2[1] + (-49./75.)*g2[0])*Yd2Tr + ((-63./5.)*g2[1] + (-81./25.)*g2[0])*Ye2Tr
                                    + (484./15.)*g4[2] + (-1096./75.)*g2[0]*g2[2] + (-24./5.)*g2[1]*g2[2] + (-69./25.)*g2[0]*g2[1] + (-81./5.)*g4[1] + (-32117./375.)*g4[0]).real();
        dX.g[1] += loopfactor3*g3[1]*((24.)*Yu4Tr + (18.)*Yu2Tr*Yu2Tr + (12.)*Yu2Yd2Tr + (24.)*Yd4Tr + (18.)*Yd2Tr*Yd2Tr + (12.)*Ye2Tr*Yd2Tr + (8.)*Ye4Tr + (2.)*Ye2Tr*Ye2Tr
                                    + ((-32.)*g2[2] + (-33.)*g2[1] + (-29./5.)*g2[0])*Yu2Tr + ((-32.)*g2[2] + (-33.)*g2[1] + (-11./5.)*g2[0])*Yd2Tr + ((-11.)*g2[1] + (-21./5.)*g2[0])*Ye2Tr
                                    + (44.)*g4[2] + (9./5.)*g2[0]*g2[1] + (24.)*g2[1]*g2[2] + (-8./5.)*g2[0]*g2[2] + (35.)*g4[1] + (-457./25.)*g4[0]).real();
        dX.g[2] += loopfactor3*g3[2]*((12.)*Yu4Tr + (18.)*Yu2Tr*Yu2Tr + (8.)*Yu2Yd2Tr + (12.)*Yd4Tr + (18.)*Yd2Tr*Yd2Tr + (6.)*Ye2Tr*Yd2Tr + ((-104./3.)*g2[2] + (-12.)*g2[1])*(Yu2Tr + Yd2Tr)
                                    + ((-44./15.)*Yu2Tr + (-32./15.)*Yd2Tr)*g2[0] + (347./3.)*g4[2] + (6.)*g2[1]*g2[2] + (22./15.)*g2[0]*g2[2] + (-3./5.)*g2[0]*g2[1] + (-27.)*g4[1] + (-1702./75.)*g4[0]).real();
      }
    }
  } else {
    dX.setZero();
  }
};
