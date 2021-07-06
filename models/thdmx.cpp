#include "thdm.h"

#define loopfactor 0.006332573977646111 // 1/(4 Pi)^2
#define loopfactor2 0.00004010149318236069 // 1/(4 Pi)^4
#define loopfactor3 0.00000025394567219137 // 1/(4 Pi)^6
#define z3 1.2020569031595942854 // Zeta(3)

using namespace Eigen;

// beta functions for 2HDM Type X
void thdmx::operator()(const thdmx &X, thdmx &dX, const double) {
  if(check()) {

    // variables for powers of parameters for nloops >= 1
    gauge<3> g2(X.g.array().square().matrix());
    gauge<3> g3(X.g.array().cube().matrix());
    gauge<3> g4(g2.array().square().matrix());
    yukawa<3,3> Yu2 = X.Yu.adjoint()*X.Yu; std::complex<double> Yu2Tr = Yu2.trace();
    yukawa<3,3> Yd2 = X.Yd.adjoint()*X.Yd; std::complex<double> Yd2Tr = Yd2.trace();
    yukawa<3,3> Ye2 = X.Ye.adjoint()*X.Ye; std::complex<double> Ye2Tr = Ye2.trace();
    yukawa<3,3> Yu4 = Yu2*Yu2;             std::complex<double> Yu4Tr = Yu4.trace();
    yukawa<3,3> Yd4 = Yd2*Yd2;             std::complex<double> Yd4Tr = Yd4.trace();
    yukawa<3,3> Ye4 = Ye2*Ye2;             std::complex<double> Ye4Tr = Ye4.trace();
    yukawa<3,3> Yu2Yd2 = Yu2*Yd2;          std::complex<double> Yu2Yd2Tr = Yu2Yd2.trace();
    yukawa<3,3> Yd2Yu2 = Yd2*Yu2;
    self<5> La2(X.La.array().square().matrix());

    dX.g[0] = loopfactor*((21./5.)*g3[0]);
    dX.g[1] = loopfactor*((-3.)*g3[1]);
    dX.g[2] = loopfactor*((-7.)*g3[2]);
    
    if(!X.weylordering || (X.weylordering && X.nloops > 1)) {
    dX.Yu = loopfactor*((-3./2.)*X.Yu*Yd2 + (3./2.)*X.Yu*Yu2 + (3.)*X.Yu*Yd2Tr + (3.)*X.Yu*Yu2Tr + (-17./20.)*X.Yu*g2[0] + (-9./4.)*X.Yu*g2[1] + (-8.)*X.Yu*g2[2]);
    dX.Yd = loopfactor*((3./2.)*X.Yd*Yd2 + (-3./2.)*X.Yd*Yu2 + (3.)*X.Yd*Yd2Tr + (3.)*X.Yd*Yu2Tr + (-1./4.)*X.Yd*g2[0] + (-9./4.)*X.Yd*g2[1] + (-8.)*X.Yd*g2[2]);
    dX.Ye = loopfactor*((3./2.)*X.Ye*Ye2 + X.Ye*Ye2Tr + (-9./4.)*X.Ye*g2[0] + (-9./4.)*X.Ye*g2[1]);
    }
    
    if(!X.weylordering || (X.weylordering && X.nloops > 2)) {
    dX.La[0] = loopfactor*((-8.)*Ye4Tr + (9./5.)*g2[0]*g2[1] + (27./50.)*g4[0] + (9./2.)*g4[1] + (4.)*Ye2Tr*X.La[0] + (-9./5.)*g2[0]*X.La[0] + (-9.)*g2[1]*X.La[0]
                         + (2.)*X.La[2]*X.La[3] + std::conj(La[4])*X.La[4] + (6.)*La2[0] + (2.)*La2[2] + La2[3]).real();
    dX.La[1] = loopfactor*((-24.)*Yd4Tr + (-24.)*Yu4Tr + (9./5.)*g2[0]*g2[1] + (27./50.)*g4[0] + (9./2.)*g4[1] + (12.)*Yd2Tr*X.La[1] + (12.)*Yu2Tr*X.La[1] + (-9./5.)*g2[0]*X.La[1]
                         + (-9.)*g2[1]*X.La[1] + (2.)*X.La[2]*X.La[3] + std::conj(La[4])*X.La[4] + (6.)*La2[1] + (2.)*La2[2] + La2[3]).real();
    dX.La[2] = loopfactor*((-9./5.)*g2[0]*g2[1] + (27./50.)*g4[0] + (9./2.)*g4[1] + (6.)*Yd2Tr*X.La[2] + (2.)*Ye2Tr*X.La[2] + (6.)*Yu2Tr*X.La[2] + (-9./5.)*g2[0]*X.La[2]
                         + (-9.)*g2[1]*X.La[2] + (3.)*X.La[0]*X.La[2] + (3.)*X.La[1]*X.La[2] + X.La[0]*X.La[3] + X.La[1]*X.La[3] + std::conj(La[4])*X.La[4] + (2.)*La2[2] + La2[3]).real();
    dX.La[3] = loopfactor*((18./5.)*g2[0]*g2[1] + (6.)*Yd2Tr*X.La[3] + (2.)*Ye2Tr*X.La[3] + (6.)*Yu2Tr*X.La[3] + (-9./5.)*g2[0]*X.La[3] + (-9.)*g2[1]*X.La[3] + X.La[0]*X.La[3] + X.La[1]*X.La[3]
                         + (4.)*X.La[2]*X.La[3] + (4.)*std::conj(La[4])*X.La[4] + (2.)*La2[3]).real();
    dX.La[4] = loopfactor*((6.)*Yd2Tr*X.La[4] + (2.)*Ye2Tr*X.La[4] + (6.)*Yu2Tr*X.La[4] + (-9./5.)*g2[0]*X.La[4] + (-9.)*g2[1]*X.La[4] + X.La[0]*X.La[4] + X.La[1]*X.La[4] + (4.)*X.La[2]*X.La[4] + (6.)*X.La[3]*X.La[4]);
    }

    if(nloops > 1) {
      gauge<3> g5(g3.cwiseProduct(g2));
      gauge<3> g6(g3.array().square().matrix());
      self<5> La3(La2.cwiseProduct(La));
      yukawa<3,3> Yu6 = Yu4*Yu2;             std::complex<double> Yu6Tr = Yu6.trace();
      yukawa<3,3> Yd6 = Yd4*Yd2;             std::complex<double> Yd6Tr = Yd6.trace();
      yukawa<3,3> Ye6 = Ye4*Ye2;             std::complex<double> Ye6Tr = Ye6.trace();
      yukawa<3,3> Yu4Yd2 = Yu4*Yd2;          std::complex<double> Yu4Yd2Tr = Yu4Yd2.trace();
      yukawa<3,3> Yu2Yd4 = Yu2*Yd4;          std::complex<double> Yu2Yd4Tr = Yu2Yd4.trace();

      dX.g[0] += loopfactor2*((-1./2.)*Yd2Tr*g3[0] + (-3./2.)*Ye2Tr*g3[0] + (-17./10.)*Yu2Tr*g3[0] + (18./5.)*g2[1]*g3[0] + (44./5.)*g2[2]*g3[0] + (104./25.)*g5[0]).real();
      dX.g[1] += loopfactor2*((-3./2.)*Yd2Tr*g3[1] + (-1./2.)*Ye2Tr*g3[1] + (-3./2.)*Yu2Tr*g3[1] + (6./5.)*g2[0]*g3[1] + (12.)*g2[2]*g3[1] + (8.)*g5[1]).real();
      dX.g[2] += loopfactor2*((-2.)*Yd2Tr*g3[2] + (-2.)*Yu2Tr*g3[2] + (11./10.)*g2[0]*g3[2] + (9./2.)*g2[1]*g3[2] + (-26.)*g5[2]).real();
      
      if(!X.weylordering || (X.weylordering && X.nloops > 2)) {
      dX.Yu += loopfactor2*((-1./4.)*X.Yu*Yd2Yu2 + (11./4.)*X.Yu*Yd4 + (-1.)*X.Yu*Yu2Yd2 + (3./2.)*X.Yu*Yu4 + (15./4.)*X.Yu*Yd2*Yd2Tr + (-27./4.)*X.Yu*Yu2*Yd2Tr + (-27./4.)*X.Yu*Yd4Tr
                          + (15./4.)*X.Yu*Yd2*Yu2Tr + (-27./4.)*X.Yu*Yu2*Yu2Tr + (3./2.)*X.Yu*Yu2Yd2Tr + (-27./4.)*X.Yu*Yu4Tr + (-43./80.)*X.Yu*Yd2*g2[0] + (223./80.)*X.Yu*Yu2*g2[0]
                          + (5./8.)*X.Yu*Yd2Tr*g2[0] + (17./8.)*X.Yu*Yu2Tr*g2[0] + (9./16.)*X.Yu*Yd2*g2[1] + (135./16.)*X.Yu*Yu2*g2[1] + (45./8.)*X.Yu*Yd2Tr*g2[1]
                          + (45./8.)*X.Yu*Yu2Tr*g2[1] + (-9./20.)*X.Yu*g2[0]*g2[1] + (-16.)*X.Yu*Yd2*g2[2] + (16.)*X.Yu*Yu2*g2[2] + (20.)*X.Yu*Yd2Tr*g2[2] + (20.)*X.Yu*Yu2Tr*g2[2]
                          + (19./15.)*X.Yu*g2[0]*g2[2] + (9.)*X.Yu*g2[1]*g2[2] + (1267./600.)*X.Yu*g4[0] +  (-21./4.)*X.Yu*g4[1] + (-108.)*X.Yu*g4[2] + (-6.)*X.Yu*Yu2*X.La[1]
                          + X.Yu*X.La[2]*X.La[3] + (3./2.)*std::conj(La[4])*X.Yu*X.La[4] + (3./2.)*X.Yu*La2[1] + X.Yu*La2[2] + X.Yu*La2[3]);
      dX.Yd += loopfactor2*((-1.)*X.Yd*Yd2Yu2 + (3./2.)*X.Yd*Yd4 + (-1./4.)*X.Yd*Yu2Yd2 + (11./4.)*X.Yd*Yu4 + (-27./4.)*X.Yd*Yd2*Yd2Tr + (15./4.)*X.Yd*Yu2*Yd2Tr + (-27./4.)*X.Yd*Yd4Tr
                          + (-27./4.)*X.Yd*Yd2*Yu2Tr + (15./4.)*X.Yd*Yu2*Yu2Tr + (3./2.)*X.Yd*Yu2Yd2Tr + (-27./4.)*X.Yd*Yu4Tr + (187./80.)*X.Yd*Yd2*g2[0] + (-79./80.)*X.Yd*Yu2*g2[0]
                          + (5./8.)*X.Yd*Yd2Tr*g2[0] + (17./8.)*X.Yd*Yu2Tr*g2[0] + (135./16.)*X.Yd*Yd2*g2[1] + (9./16.)*X.Yd*Yu2*g2[1] + (45./8.)*X.Yd*Yd2Tr*g2[1] + (45./8.)*X.Yd*Yu2Tr*g2[1]
                          + (-27./20.)*X.Yd*g2[0]*g2[1] + (16.)*X.Yd*Yd2*g2[2] + (-16.)*X.Yd*Yu2*g2[2] + (20.)*X.Yd*Yd2Tr*g2[2] + (20.)*X.Yd*Yu2Tr*g2[2] + (31./15.)*X.Yd*g2[0]*g2[2]
                          + (9.)*X.Yd*g2[1]*g2[2] + (-113./600.)*X.Yd*g4[0] + (-21./4.)*X.Yd*g4[1] + (-108.)*X.Yd*g4[2] + (-6.)*X.Yd*Yd2*X.La[1] + X.Yd*X.La[2]*X.La[3]
                          + (3./2.)*std::conj(La[4])*X.Yd*X.La[4] + (3./2.)*X.Yd*La2[1] + X.Yd*La2[2] + X.Yd*La2[3]);
      dX.Ye += loopfactor2*((3./2.)*X.Ye*Ye4 + (-9./4.)*X.Ye*Ye2*Ye2Tr + (-9./4.)*X.Ye*Ye4Tr + (387./80.)*X.Ye*Ye2*g2[0] + (15./8.)*X.Ye*Ye2Tr*g2[0] + (135./16.)*X.Ye*Ye2*g2[1]
                          + (15./8.)*X.Ye*Ye2Tr*g2[1] + (27./20.)*X.Ye*g2[0]*g2[1] + (1449./200.)*X.Ye*g4[0] + (-21./4.)*X.Ye*g4[1] + (-6.)*X.Ye*Ye2*X.La[0] + X.Ye*X.La[2]*X.La[3]
                          + (3./2.)*std::conj(La[4])*X.Ye*X.La[4] + (3./2.)*X.Ye*La2[0] + X.Ye*La2[2] + X.Ye*La2[3]);
      }
      
      if(!X.weylordering || (X.weylordering && X.nloops > 3)) {
      dX.La[0] += loopfactor2*((-5.)*X.La[0]*X.La[2]*X.La[3] + (-7./2.)*std::conj(X.La[4])*X.La[0]*X.La[4] + (-10.)*std::conj(X.La[4])*X.La[2]*X.La[4] + (-11.)*std::conj(X.La[4])*X.La[3]*X.La[4]
                              + (12./5.)*X.La[2]*X.La[3]*g2[0] + (-3./5.)*std::conj(X.La[4])*X.La[4]*g2[0] + (12.)*X.La[2]*X.La[3]*g2[1] + (117./20.)*X.La[0]*g2[0]*g2[1] + (3.)*X.La[3]*g2[0]*g2[1]
                              + (1953./200.)*X.La[0]*g4[0] + (9./5.)*X.La[2]*g4[0] + (9./10.)*X.La[3]*g4[0] + (-1719./100.)*g2[1]*g4[0] + (-51./8.)*X.La[0]*g4[1] + (15.)*X.La[2]*g4[1] + (15./2.)*X.La[3]*g4[1]
                              + (-303./20.)*g2[0]*g4[1] + (-3537./500.)*g6[0] + (291./4.)*g6[1] + (27./5.)*g2[0]*La2[0] + (27.)*g2[1]*La2[0] + (-5.)*X.La[0]*La2[2] + (-6.)*X.La[3]*La2[2] + (12./5.)*g2[0]*La2[2]
                              + (12.)*g2[1]*La2[2] + (-3.)*X.La[0]*La2[3] + (-8.)*X.La[2]*La2[3] + (6./5.)*g2[0]*La2[3] + (3.)*g2[1]*La2[3] + (-39./2.)*La3[0] + (-4.)*La3[2] + (-3.)*La3[3]
                              + (-12.)*X.La[2]*X.La[3]*Yd2Tr + (-6.)*std::conj(X.La[4])*X.La[4]*Yd2Tr + (-1.)*X.La[0]*Ye4Tr + (40.)*Ye6Tr + (-12.)*X.La[2]*X.La[3]*Yu2Tr + (-6.)*std::conj(X.La[4])*X.La[4]*Yu2Tr
                              + (15./2.)*X.La[0]*Ye2Tr*g2[0] + (-48./5.)*Ye4Tr*g2[0] + (15./2.)*X.La[0]*Ye2Tr*g2[1] + (66./5.)*Ye2Tr*g2[0]*g2[1] + (-9.)*Ye2Tr*g4[0] + (-3.)*Ye2Tr*g4[1] + (-12.)*Ye2Tr*La2[0]
                              + (-12.)*Yd2Tr*La2[2] + (-12.)*Yu2Tr*La2[2] + (-6.)*Yd2Tr*La2[3] + (-6.)*Yu2Tr*La2[3]).real();
      dX.La[1] += loopfactor2*((-5.)*X.La[1]*X.La[2]*X.La[3] + (-7./2.)*std::conj(X.La[4])*X.La[1]*X.La[4] + (-10.)*std::conj(X.La[4])*X.La[2]*X.La[4] + (-11.)*std::conj(X.La[4])*X.La[3]*X.La[4]
                              + (12./5.)*X.La[2]*X.La[3]*g2[0] + (-3./5.)*std::conj(X.La[4])*X.La[4]*g2[0] + (12.)*X.La[2]*X.La[3]*g2[1] + (117./20.)*X.La[1]*g2[0]*g2[1] + (3.)*X.La[3]*g2[0]*g2[1]
                              + (1953./200.)*X.La[1]*g4[0] + (9./5.)*X.La[2]*g4[0] + (9./10.)*X.La[3]*g4[0] + (-1719./100.)*g2[1]*g4[0] + (-51./8.)*X.La[1]*g4[1] + (15.)*X.La[2]*g4[1] + (15./2.)*X.La[3]*g4[1]
                              + (-303./20.)*g2[0]*g4[1] + (-3537./500.)*g6[0] + (291./4.)*g6[1] + (27./5.)*g2[0]*La2[1] + (27.)*g2[1]*La2[1] + (-5.)*X.La[1]*La2[2] + (-6.)*X.La[3]*La2[2] + (12./5.)*g2[0]*La2[2]
                              + (12.)*g2[1]*La2[2] + (-3.)*X.La[1]*La2[3] + (-8.)*X.La[2]*La2[3] + (6./5.)*g2[0]*La2[3] + (3.)*g2[1]*La2[3] + (-39./2.)*La3[1] + (-4.)*La3[2] + (-3.)*La3[3]
                              + (-3.)*X.La[1]*Yd4Tr + (120.)*Yd6Tr + (-4.)*X.La[2]*X.La[3]*Ye2Tr + (-2.)*std::conj(X.La[4])*X.La[4]*Ye2Tr + (-42.)*X.La[1]*Yu2Yd2Tr + (-24.)*Yu2Yd4Tr + (-3.)*X.La[1]*Yu4Tr
                              + (-24.)*Yu4Yd2Tr + (120.)*Yu6Tr + (5./2.)*X.La[1]*Yd2Tr*g2[0] + (16./5.)*Yd4Tr*g2[0] + (17./2.)*X.La[1]*Yu2Tr*g2[0] + (-32./5.)*Yu4Tr*g2[0] + (45./2.)*X.La[1]*Yd2Tr*g2[1]
                              + (45./2.)*X.La[1]*Yu2Tr*g2[1] + (54./5.)*Yd2Tr*g2[0]*g2[1] + (126./5.)*Yu2Tr*g2[0]*g2[1] + (80.)*X.La[1]*Yd2Tr*g2[2] + (-128.)*Yd4Tr*g2[2] + (80.)*X.La[1]*Yu2Tr*g2[2] + (-128.)*Yu4Tr*g2[2]
                              + (9./5.)*Yd2Tr*g4[0] + (-171./25.)*Yu2Tr*g4[0] + (-9.)*Yd2Tr*g4[1] + (-9.)*Yu2Tr*g4[1] + (-36.)*Yd2Tr*La2[1] + (-36.)*Yu2Tr*La2[1] + (-4.)*Ye2Tr*La2[2] + (-2.)*Ye2Tr*La2[3]).real();
      dX.La[2] += loopfactor2*((-4.)*X.La[0]*X.La[2]*X.La[3] + (-4.)*X.La[1]*X.La[2]*X.La[3] + (-9./2.)*std::conj(X.La[4])*X.La[0]*X.La[4] + (-9./2.)*std::conj(X.La[4])*X.La[1]*X.La[4]
                              + (-9./2.)*std::conj(X.La[4])*X.La[2]*X.La[4] + (-11.)*std::conj(X.La[4])*X.La[3]*X.La[4] + (18./5.)*X.La[0]*X.La[2]*g2[0] + (18./5.)*X.La[1]*X.La[2]*g2[0]
                              + (6./5.)*X.La[0]*X.La[3]*g2[0] + (6./5.)*X.La[1]*X.La[3]*g2[0] + (6./5.)*std::conj(X.La[4])*X.La[4]*g2[0] + (18.)*X.La[0]*X.La[2]*g2[1] + (18.)*X.La[1]*X.La[2]*g2[1]
                              + (9.)*X.La[0]*X.La[3]*g2[1] + (9.)*X.La[1]*X.La[3]*g2[1] + (-6.)*X.La[2]*X.La[3]*g2[1] + (-3./2.)*X.La[0]*g2[0]*g2[1] + (-3./2.)*X.La[1]*g2[0]*g2[1] + (33./20.)*X.La[2]*g2[0]*g2[1]
                              + (-9./5.)*X.La[3]*g2[0]*g2[1] + (27./20.)*X.La[0]*g4[0] + (27./20.)*X.La[1]*g4[0] + (1773./200.)*X.La[2]*g4[0] + (9./10.)*X.La[3]*g4[0] + (909./100.)*g2[1]*g4[0] + (45./4.)*X.La[0]*g4[1]
                              + (45./4.)*X.La[1]*g4[1] + (-111./8.)*X.La[2]*g4[1] + (15./2.)*X.La[3]*g4[1] + (33./20.)*g2[0]*g4[1] + (-3537./500.)*g6[0] + (291./4.)*g6[1] + (-15./4.)*X.La[2]*La2[0]
                              + (-1.)*X.La[3]*La2[0] + (-15./4.)*X.La[2]*La2[1] + (-1.)*X.La[3]*La2[1] + (-9.)*X.La[0]*La2[2] + (-9.)*X.La[1]*La2[2] + (-1.)*X.La[3]*La2[2] + (3./5.)*g2[0]*La2[2] + (3.)*g2[1]*La2[2]
                              + (-7./2.)*X.La[0]*La2[3] + (-7./2.)*X.La[1]*La2[3] + (-4.)*X.La[2]*La2[3] + (-3./5.)*g2[0]*La2[3] + (3.)*g2[1]*La2[3] + (-3.)*La3[2] + (-3.)*La3[3]
                              + (-18.)*X.La[1]*X.La[2]*Yd2Tr + (-6.)*X.La[1]*X.La[3]*Yd2Tr + (-3.)*std::conj(X.La[4])*X.La[4]*Yd2Tr + (-27./2.)*X.La[2]*Yd4Tr + (-6.)*X.La[1]*X.La[2]*Ye2Tr
                              + (-2.)*X.La[1]*X.La[3]*Ye2Tr + (-1.)*std::conj(X.La[4])*X.La[4]*Ye2Tr + (-9./2.)*X.La[2]*Ye4Tr + (-21.)*X.La[2]*Yu2Yd2Tr + (-24.)*X.La[3]*Yu2Yd2Tr + (-18.)*X.La[1]*X.La[2]*Yu2Tr
                              + (-6.)*X.La[1]*X.La[3]*Yu2Tr + (-3.)*std::conj(X.La[4])*X.La[4]*Yu2Tr + (-27./2.)*X.La[2]*Yu4Tr + (5./4.)*X.La[2]*Yd2Tr*g2[0] + (15./4.)*X.La[2]*Ye2Tr*g2[0]
                              + (17./4.)*X.La[2]*Yu2Tr*g2[0] + (45./4.)*X.La[2]*Yd2Tr*g2[1] + (15./4.)*X.La[2]*Ye2Tr*g2[1] + (45./4.)*X.La[2]*Yu2Tr*g2[1] + (-27./5.)*Yd2Tr*g2[0]*g2[1] + (-33./5.)*Ye2Tr*g2[0]*g2[1]
                              + (-63./5.)*Yu2Tr*g2[0]*g2[1] + (40.)*X.La[2]*Yd2Tr*g2[2] + (40.)*X.La[2]*Yu2Tr*g2[2] + (9./10.)*Yd2Tr*g4[0] + (-9./2.)*Ye2Tr*g4[0] + (-171./50.)*Yu2Tr*g4[0] + (-9./2.)*Yd2Tr*g4[1]
                              + (-3./2.)*Ye2Tr*g4[1] + (-9./2.)*Yu2Tr*g4[1] + (-6.)*Yd2Tr*La2[2] + (-2.)*Ye2Tr*La2[2] + (-6.)*Yu2Tr*La2[2] + (-3.)*Yd2Tr*La2[3] + (-1.)*Ye2Tr*La2[3] + (-3.)*Yu2Tr*La2[3]
                              + (6.)*X.La[0]*X.La[2]*Ye2Tr + (-6.)*X.La[1]*X.La[2]*Ye2Tr + (2.)*X.La[0]*X.La[3]*Ye2Tr + (-2.)*X.La[1]*X.La[3]*Ye2Tr).real();
      dX.La[3] += loopfactor2*((-10.)*X.La[0]*X.La[2]*X.La[3] + (-10.)*X.La[1]*X.La[2]*X.La[3] + (-6.)*std::conj(X.La[4])*X.La[0]*X.La[4] + (-6.)*std::conj(X.La[4])*X.La[1]*X.La[4]
                              + (-12.)*std::conj(X.La[4])*X.La[2]*X.La[4] + (-13./2.)*std::conj(X.La[4])*X.La[3]*X.La[4] + (6./5.)*X.La[0]*X.La[3]*g2[0] + (6./5.)*X.La[1]*X.La[3]*g2[0]
                              + (6./5.)*X.La[2]*X.La[3]*g2[0] + (24./5.)*std::conj(X.La[4])*X.La[4]*g2[0] + (18.)*X.La[2]*X.La[3]*g2[1] + (27.)*std::conj(X.La[4])*X.La[4]*g2[1] + (3.)*X.La[0]*g2[0]*g2[1]
                              + (3.)*X.La[1]*g2[0]*g2[1] + (6./5.)*X.La[2]*g2[0]*g2[1] + (153./20.)*X.La[3]*g2[0]*g2[1] + (1413./200.)*X.La[3]*g4[0] + (-657./25.)*g2[1]*g4[0] + (-231./8.)*X.La[3]*g4[1]
                              + (-84./5.)*g2[0]*g4[1] + (-7./4.)*X.La[3]*La2[0] + (-7./4.)*X.La[3]*La2[1] + (-7.)*X.La[3]*La2[2] + (-5.)*X.La[0]*La2[3] + (-5.)*X.La[1]*La2[3] + (-7.)*X.La[2]*La2[3]
                              + (12./5.)*g2[0]*La2[3] + (9.)*g2[1]*La2[3]
                              + (-6.)*X.La[1]*X.La[3]*Yd2Tr + (-12.)*X.La[2]*X.La[3]*Yd2Tr + (-12.)*std::conj(X.La[4])*X.La[4]*Yd2Tr + (-27./2.)*X.La[3]*Yd4Tr + (-2.)*X.La[1]*X.La[3]*Ye2Tr
                              + (-4.)*X.La[2]*X.La[3]*Ye2Tr + (-4.)*std::conj(X.La[4])*X.La[4]*Ye2Tr + (-9./2.)*X.La[3]*Ye4Tr + (-6.)*X.La[1]*X.La[3]*Yu2Tr + (-12.)*X.La[2]*X.La[3]*Yu2Tr
                              + (-12.)*std::conj(X.La[4])*X.La[4]*Yu2Tr + (27.)*X.La[3]*Yu2Yd2Tr + (-27./2.)*X.La[3]*Yu4Tr + (5./4.)*X.La[3]*Yd2Tr*g2[0] + (15./4.)*X.La[3]*Ye2Tr*g2[0]
                              + (17./4.)*X.La[3]*Yu2Tr*g2[0] + (45./4.)*X.La[3]*Yd2Tr*g2[1] + (15./4.)*X.La[3]*Ye2Tr*g2[1] + (45./4.)*X.La[3]*Yu2Tr*g2[1] + (54./5.)*Yd2Tr*g2[0]*g2[1] + (66./5.)*Ye2Tr*g2[0]*g2[1]
                              + (126./5.)*Yu2Tr*g2[0]*g2[1] + (40.)*X.La[3]*Yd2Tr*g2[2] + (40.)*X.La[3]*Yu2Tr*g2[2] + (-6.)*Yd2Tr*La2[3] + (-2.)*Ye2Tr*La2[3] + (-6.)*Yu2Tr*La2[3]
                              + (-2.)*X.La[0]*X.La[3]*Ye2Tr + (2.)*X.La[1]*X.La[3]*Ye2Tr).real();
      dX.La[4] += loopfactor2*((-10.)*X.La[0]*X.La[2]*X.La[4] + (-10.)*X.La[1]*X.La[2]*X.La[4] + (-11.)*X.La[0]*X.La[3]*X.La[4] + (-11.)*X.La[1]*X.La[3]*X.La[4] + (-19.)*X.La[2]*X.La[3]*X.La[4]
                              + (-3./5.)*X.La[0]*X.La[4]*g2[0] + (-3./5.)*X.La[1]*X.La[4]*g2[0] + (24./5.)*X.La[2]*X.La[4]*g2[0] + (36./5.)*X.La[3]*X.La[4]*g2[0] + (18.)*X.La[2]*X.La[4]*g2[1]
                              + (36.)*X.La[3]*X.La[4]*g2[1] + (57./20.)*X.La[4]*g2[0]*g2[1] + (1413./200.)*X.La[4]*g4[0] + (-231./8.)*X.La[4]*g4[1] + (-7./4.)*X.La[4]*La2[0] + (-7./4.)*X.La[4]*La2[1]
                              + (-7.)*X.La[4]*La2[2] + (-8.)*X.La[4]*La2[3] + (3./2.)*std::conj(X.La[4])*La2[4]
                              + (-6.)*X.La[1]*X.La[4]*Yd2Tr + (-12.)*X.La[2]*X.La[4]*Yd2Tr + (-18.)*X.La[3]*X.La[4]*Yd2Tr + (-3./2.)*X.La[4]*Yd4Tr + (-2.)*X.La[1]*X.La[4]*Ye2Tr + (-4.)*X.La[2]*X.La[4]*Ye2Tr
                              + (-6.)*X.La[3]*X.La[4]*Ye2Tr + (-1./2.)*X.La[4]*Ye4Tr + (3.)*X.La[4]*Yu2Yd2Tr + (-6.)*X.La[1]*X.La[4]*Yu2Tr + (-12.)*X.La[2]*X.La[4]*Yu2Tr + (-18.)*X.La[3]*X.La[4]*Yu2Tr
                              + (-3./2.)*X.La[4]*Yu4Tr + (5./4.)*X.La[4]*Yd2Tr*g2[0] + (15./4.)*X.La[4]*Ye2Tr*g2[0] + (17./4.)*X.La[4]*Yu2Tr*g2[0] + (45./4.)*X.La[4]*Yd2Tr*g2[1] + (15./4.)*X.La[4]*Ye2Tr*g2[1]
                              + (45./4.)*X.La[4]*Yu2Tr*g2[1] + (40.)*X.La[4]*Yd2Tr*g2[2] + (40.)*X.La[4]*Yu2Tr*g2[2] + (-2.)*X.La[0]*X.La[4]*Ye2Tr + (2.)*X.La[1]*X.La[4]*Ye2Tr);
      }

      if(nloops > 2) {
        gauge<3> g7(g5.cwiseProduct(g2));
        gauge<3> g8(g4.array().square().matrix());
        std::complex<double> Yu2Tr2 = Yu2Tr*Yu2Tr;
        std::complex<double> Yd2Tr2 = Yd2Tr*Yd2Tr;
        std::complex<double> Ye2Tr2 = Ye2Tr*Ye2Tr;
        std::complex<double> Yu4Tr2 = Yu4Tr*Yu4Tr;
        std::complex<double> Yd4Tr2 = Yd4Tr*Yd4Tr;
        std::complex<double> Ye4Tr2 = Ye4Tr*Ye4Tr;
        std::complex<double> Yu2Yd2Tr2 = Yu2Yd2Tr*Yu2Yd2Tr;
        yukawa<3,3> Yd2Yu4 = Yd2*Yu4;                yukawa<3,3> Yu2Yd2Yu2 = Yu2*Yd2*Yu2;
        yukawa<3,3> Yd4Yu2 = Yd4*Yu2;                yukawa<3,3> Yd2Yu2Yd2 = Yd2*Yu2*Yd2;
        yukawa<3,3> Yu8 = Yu6*Yu2;                   std::complex<double> Yu8Tr = Yu8.trace();
        yukawa<3,3> Yd8 = Yd6*Yd2;                   std::complex<double> Yd8Tr = Yd8.trace();
        yukawa<3,3> Ye8 = Ye6*Ye2;                   std::complex<double> Ye8Tr = Ye8.trace();
        yukawa<3,3> Yu4Yd4 = Yu4*Yd4;                std::complex<double> Yu4Yd4Tr = Yu4Yd4.trace();
        yukawa<3,3> Yu2Yd6 = Yu2*Yd6;                std::complex<double> Yu2Yd6Tr = Yu2Yd6.trace();
        yukawa<3,3> Yu6Yd2 = Yu6*Yd2;                std::complex<double> Yu6Yd2Tr = Yu6Yd2.trace();
        yukawa<3,3> Yu2Yd2Yu2Yd2 = Yu2*Yd2*Yu2*Yd2;  std::complex<double> Yu2Yd2Yu2Yd2Tr = Yu2Yd2Yu2Yd2.trace();

        dX.g[0] += loopfactor3*((-6./5.)*X.La[2]*X.La[3]*g3[0] + (-9./10.)*std::conj(X.La[4])*X.La[4]*g3[0] + (51./40.)*Yd2Tr2*g3[0] + (183./80.)*Yd4Tr*g3[0] + (99./40.)*Ye2Tr2*g3[0] + (261./80.)*Ye4Tr*g3[0]
                              + (177./20.)*Yd2Tr*Yu2Tr*g3[0] + (303./40.)*Yu2Tr2*g3[0] + (3./8.)*Yu2Yd2Tr*g3[0] + (339./80.)*Yu4Tr*g3[0] + (9./20.)*X.La[0]*g2[1]*g3[0] + (9./20.)*X.La[1]*g2[1]*g3[0]
                              + (9./10.)*X.La[3]*g2[1]*g3[0] + (-1311./160.)*Yd2Tr*g2[1]*g3[0] + (-1629./160.)*Ye2Tr*g2[1]*g3[0] + (-471./32.)*Yu2Tr*g2[1]*g3[0] + (-17./5.)*Yd2Tr*g2[2]*g3[0] + (-29./5.)*Yu2Tr*g2[2]*g3[0]
                              + (-3./5.)*g2[1]*g2[2]*g3[0] + (2883./160.)*g3[0]*g4[1] + (297./5.)*g3[0]*g4[2] + (27./100.)*X.La[0]*g5[0] + (27./100.)*X.La[1]*g5[0] + (9./25.)*X.La[2]*g5[0] + (9./50.)*X.La[3]*g5[0]
                              + (-1267./800.)*Yd2Tr*g5[0] + (-2529./800.)*Ye2Tr*g5[0] + (-2827./800.)*Yu2Tr*g5[0] + (699./400.)*g2[1]*g5[0] + (-137./75.)*g2[2]*g5[0] + (-42439./2400.)*g7[0] + (-9./20.)*g3[0]*La2[0]
                              + (-9./20.)*g3[0]*La2[1] + (-3./10.)*g3[0]*La2[2] + (-3./10.)*g3[0]*La2[3]).real();
        dX.g[1] += loopfactor3*((-3./2.)*std::conj(X.La[4])*X.La[4]*g3[1] + (45./8.)*Yd2Tr2*g3[1] + (57./16.)*Yd4Tr*g3[1] + (5./8.)*Ye2Tr2*g3[1] + (19./16.)*Ye4Tr*g3[1] + (45./4.)*Yd2Tr*Yu2Tr*g3[1] + (45./8.)*Yu2Tr2*g3[1]
                              + (27./8.)*Yu2Yd2Tr*g3[1] + (57./16.)*Yu4Tr*g3[1] + (3./20.)*X.La[0]*g2[0]*g3[1] + (3./20.)*X.La[1]*g2[0]*g3[1] + (3./10.)*X.La[3]*g2[0]*g3[1] + (-533./160.)*Yd2Tr*g2[0]*g3[1]
                              + (-51./32.)*Ye2Tr*g2[0]*g3[1] + (-593./160.)*Yu2Tr*g2[0]*g3[1] + (-7.)*Yd2Tr*g2[2]*g3[1] + (-7.)*Yu2Tr*g2[2]*g3[1] + (-1./5.)*g2[0]*g2[2]*g3[1] + (-3907./800.)*g3[1]*g4[0] + (81.)*g3[1]*g4[2]
                              + (3./4.)*X.La[0]*g5[1] + (3./4.)*X.La[1]*g5[1] + X.La[2]*g5[1] + (1./2.)*X.La[3]*g5[1] + (-729./32.)*Yd2Tr*g5[1] + (-243./32.)*Ye2Tr*g5[1] + (-729./32.)*Yu2Tr*g5[1] + (717./80.)*g2[0]*g5[1]
                              + (39.)*g2[2]*g5[1] + (20035./96.)*g7[1] + (-2.)*g3[1]*La[2]*La[3] + (-3./4.)*g3[1]*La2[0] + (-3./4.)*g3[1]*La2[1] + (-1./2.)*g3[1]*La2[2] + (-1./2.)*g3[1]*La2[3]).real();
        dX.g[2] += loopfactor3*((21./2.)*Yd2Tr2*g3[2] + (9./2.)*Yd4Tr*g3[2] + (21.)*Yd2Tr*Yu2Tr*g3[2] + (21./2.)*Yu2Tr2*g3[2] + (-3.)*Yu2Yd2Tr*g3[2] + (9./2.)*Yu4Tr*g3[2] + (-89./40.)*Yd2Tr*g2[0]*g3[2]
                              + (-101./40.)*Yu2Tr*g2[0]*g3[2] + (-93./8.)*Yd2Tr*g2[1]*g3[2] + (-93./8.)*Yu2Tr*g2[1]*g3[2] + (-3./40.)*g2[0]*g2[1]*g3[2] + (-5483./1200.)*g3[2]*g4[0] + (195./16.)*g3[2]*g4[1] + (-40.)*Yd2Tr*g5[2]
                              + (-40.)*Yu2Tr*g5[2] + (77./15.)*g2[0]*g5[2] + (21.)*g2[1]*g5[2] + (65./2.)*g7[2]).real();
                              
        if(!X.weylordering || (X.weylordering && X.nloops > 3)) {
        dX.Yu += loopfactor3*((-3./2.)*X.La[0]*X.La[2]*X.La[3]*X.Yu + (-3.)*X.La[1]*X.La[2]*X.La[3]*X.Yu + (-3./4.)*std::conj(X.La[4])*X.La[0]*X.La[4]*X.Yu
                              + (-3./2.)*std::conj(X.La[4])*X.La[1]*X.La[4]*X.Yu + (-9./2.)*std::conj(X.La[4])*X.La[2]*X.La[4]*X.Yu + (-27./4.)*std::conj(X.La[4])*X.La[3]*X.La[4]*X.Yu
                              + (-7./8.)*X.La[2]*X.La[3]*X.Yu*Yd2 + (75./16.)*std::conj(X.La[4])*X.La[4]*X.Yu*Yd2 + (3./2.)*X.La[1]*X.Yu*Yd2Yu2 + (-37./8.)*X.Yu*Yd2Yu2Yd2 + (75./16.)*X.Yu*Yd2Yu4
                              + (-15.)*X.La[1]*X.Yu*Yd4 + (-95./8.)*X.Yu*Yd4Yu2 + (9./4.)*X.Yu*Yd6 + (15./8.)*X.La[2]*X.La[3]*X.Yu*Yu2 + (45./16.)*std::conj(X.La[4])*X.La[4]*X.Yu*Yu2 + (43./8.)*X.Yu*Yu2Yd2Yu2
                              + (-183./16.)*X.Yu*Yu2Yd4 + (63./2.)*X.La[1]*X.Yu*Yu4 + (83./16.)*X.Yu*Yu4Yd2 + (-345./16.)*X.Yu*Yu6 + (-15./4.)*X.La[2]*X.La[3]*X.Yu*Yd2Tr
                              + (-45./8.)*std::conj(X.La[4])*X.La[4]*X.Yu*Yd2Tr + (21./4.)*X.Yu*Yd2Yu2*Yd2Tr + (-69./4.)*X.Yu*Yd4*Yd2Tr + (45.)*X.La[1]*X.Yu*Yu2*Yd2Tr + (3./4.)*X.Yu*Yu2Yd2*Yd2Tr
                              + (-9./4.)*X.Yu*Yu4*Yd2Tr + (117./8.)*X.Yu*Yd2*Yd2Tr2 + (-81./8.)*X.Yu*Yu2*Yd2Tr2 + (45./2.)*X.La[1]*X.Yu*Yd4Tr + (-153./8.)*X.Yu*Yd2*Yd4Tr + (27.)*X.Yu*Yu2*Yd4Tr
                              + (54.)*X.Yu*Yd2Tr*Yd4Tr + (-75./16.)*X.Yu*Yd6Tr + (-5./2.)*X.La[2]*X.La[3]*X.Yu*Ye2Tr + (-15./4.)*std::conj(X.La[4])*X.La[4]*X.Yu*Ye2Tr + (-15./4.)*X.La[2]*X.La[3]*X.Yu*Yu2Tr
                              + (-45./8.)*std::conj(X.La[4])*X.La[4]*X.Yu*Yu2Tr + (21./4.)*X.Yu*Yd2Yu2*Yu2Tr + (-69./4.)*X.Yu*Yd4*Yu2Tr + (45.)*X.La[1]*X.Yu*Yu2*Yu2Tr + (3./4.)*X.Yu*Yu2Yd2*Yu2Tr
                              + (-9./4.)*X.Yu*Yu4*Yu2Tr + (117./4.)*X.Yu*Yd2*Yd2Tr*Yu2Tr + (-81./4.)*X.Yu*Yu2*Yd2Tr*Yu2Tr + (54.)*X.Yu*Yd4Tr*Yu2Tr + (117./8.)*X.Yu*Yd2*Yu2Tr2 + (-81./8.)*X.Yu*Yu2*Yu2Tr2
                              + (177./4.)*X.Yu*Yd2*Yu2Yd2Tr + (-9./2.)*X.Yu*Yd2Tr*Yu2Yd2Tr + (-9./2.)*X.Yu*Yu2Tr*Yu2Yd2Tr + (39./16.)*X.Yu*Yu2Yd4Tr + (45./2.)*X.La[1]*X.Yu*Yu4Tr + (-153./8.)*X.Yu*Yd2*Yu4Tr
                              + (27.)*X.Yu*Yu2*Yu4Tr + (54.)*X.Yu*Yd2Tr*Yu4Tr + (54.)*X.Yu*Yu2Tr*Yu4Tr + (39./16.)*X.Yu*Yu4Yd2Tr + (-75./16.)*X.Yu*Yu6Tr + (12.)*X.Yu*Yd4Yu2*z3 + (-9./2.)*X.Yu*Yd6*z3
                              + (12.)*X.Yu*Yu2Yd4*z3 + (9./2.)*X.Yu*Yu6*z3 + (9.)*X.Yu*Yd6Tr*z3 + (-72.)*X.Yu*Yd2*Yu2Yd2Tr*z3 + (9.)*X.Yu*Yu6Tr*z3 + (3./2.)*X.La[2]*X.La[3]*X.Yu*g2[0]
                              + (9./4.)*std::conj(X.La[4])*X.La[4]*X.Yu*g2[0] + (-199./160.)*X.Yu*Yd2Yu2*g2[0] + (1047./160.)*X.Yu*Yd4*g2[0] + (-127./20.)*X.La[1]*X.Yu*Yu2*g2[0] + (3./5.)*X.Yu*Yu2Yd2*g2[0]
                              + (-28./5.)*X.Yu*Yu4*g2[0] + (23./8.)*X.Yu*Yd2*Yd2Tr*g2[0] + (-57./10.)*X.Yu*Yu2*Yd2Tr*g2[0] + (-81./20.)*X.Yu*Yd2Tr2*g2[0] + (-909./80.)*X.Yu*Yd4Tr*g2[0] + (65./8.)*X.Yu*Yd2*Yu2Tr*g2[0]
                              + (-129./10.)*X.Yu*Yu2*Yu2Tr*g2[0] + (-81./10.)*X.Yu*Yd2Tr*Yu2Tr*g2[0] + (-81./20.)*X.Yu*Yu2Tr2*g2[0] + (-93./40.)*X.Yu*Yu2Yd2Tr*g2[0] + (-633./80.)*X.Yu*Yu4Tr*g2[0]
                              + (41./20.)*X.Yu*Yd2Yu2*z3*g2[0] + (-53./10.)*X.Yu*Yd4*z3*g2[0] + (17./20.)*X.Yu*Yu2Yd2*z3*g2[0] + (9./5.)*X.Yu*Yd2*Yd2Tr*z3*g2[0] + (-18./5.)*X.Yu*Yu2*Yd2Tr*z3*g2[0]
                              + (27./5.)*X.Yu*Yd4Tr*z3*g2[0] + (-18./5.)*X.Yu*Yd2*Yu2Tr*z3*g2[0] + (9./5.)*X.Yu*Yu2*Yu2Tr*z3*g2[0] + (24./5.)*X.Yu*Yu2Yd2Tr*z3*g2[0] + (-9./5.)*X.Yu*Yu4Tr*z3*g2[0]
                              + (15./2.)*X.La[2]*X.La[3]*X.Yu*g2[1] + (45./4.)*std::conj(X.La[4])*X.La[4]*X.Yu*g2[1] + (39./32.)*X.Yu*Yd2Yu2*g2[1] + (579./32.)*X.Yu*Yd4*g2[1] + (-135./4.)*X.La[1]*X.Yu*Yu2*g2[1]
                              + (195./16.)*X.Yu*Yu2Yd2*g2[1] + (-27./4.)*X.Yu*Yu4*g2[1] + (-135./8.)*X.Yu*Yd2*Yd2Tr*g2[1] + (-81./4.)*X.Yu*Yu2*Yd2Tr*g2[1] + (-81./4.)*X.Yu*Yd2Tr2*g2[1] + (-837./16.)*X.Yu*Yd4Tr*g2[1]
                              + (-135./8.)*X.Yu*Yd2*Yu2Tr*g2[1] + (-81./4.)*X.Yu*Yu2*Yu2Tr*g2[1] + (-81./2.)*X.Yu*Yd2Tr*Yu2Tr*g2[1] + (-81./4.)*X.Yu*Yu2Tr2*g2[1] + (-63./8.)*X.Yu*Yu2Yd2Tr*g2[1]
                              + (-837./16.)*X.Yu*Yu4Tr*g2[1] + (-9./4.)*X.Yu*Yd2Yu2*z3*g2[1] + (-45./2.)*X.Yu*Yd4*z3*g2[1] + (-9./4.)*X.Yu*Yu2Yd2*z3*g2[1] + (27.)*X.Yu*Yd2*Yd2Tr*z3*g2[1] + (-27.)*X.Yu*Yu2*Yd2Tr*z3*g2[1]
                              + (27.)*X.Yu*Yd4Tr*z3*g2[1] + (27.)*X.Yu*Yd2*Yu2Tr*z3*g2[1] + (-27.)*X.Yu*Yu2*Yu2Tr*z3*g2[1] + (27.)*X.Yu*Yu4Tr*z3*g2[1] + (117./80.)*X.La[1]*X.Yu*g2[0]*g2[1]
                              + (117./80.)*X.La[3]*X.Yu*g2[0]*g2[1] + (-1443./640.)*X.Yu*Yd2*g2[0]*g2[1] + (-2193./640.)*X.Yu*Yu2*g2[0]*g2[1] + (2589./320.)*X.Yu*Yd2Tr*g2[0]*g2[1] + (1029./64.)*X.Yu*Yu2Tr*g2[0]*g2[1]
                              + (27./10.)*X.Yu*Yd2*z3*g2[0]*g2[1] + (207./20.)*X.Yu*Yu2*z3*g2[0]*g2[1] + (81./10.)*X.Yu*Yu2Tr*z3*g2[0]*g2[1] + (28.)*X.Yu*Yd2Yu2*g2[2] + (26.)*X.Yu*Yd4*g2[2]
                              + (8.)*X.La[1]*X.Yu*Yu2*g2[2] + (-18.)*X.Yu*Yu2Yd2*g2[2] + (-76.)*X.Yu*Yu4*g2[2] + (97./2.)*X.Yu*Yd2*Yd2Tr*g2[2] + (-177./2.)*X.Yu*Yu2*Yd2Tr*g2[2] + (15./2.)*X.Yu*Yd4Tr*g2[2]
                              + (97./2.)*X.Yu*Yd2*Yu2Tr*g2[2] + (-177./2.)*X.Yu*Yu2*Yu2Tr*g2[2] + (57.)*X.Yu*Yu2Yd2Tr*g2[2] + (15./2.)*X.Yu*Yu4Tr*g2[2] + (8.)*X.Yu*Yd2Yu2*z3*g2[2] + (80.)*X.Yu*Yd4*z3*g2[2]
                              + (8.)*X.Yu*Yu2Yd2*z3*g2[2] + (-72.)*X.Yu*Yd2*Yd2Tr*z3*g2[2] + (72.)*X.Yu*Yu2*Yd2Tr*z3*g2[2] + (-72.)*X.Yu*Yd4Tr*z3*g2[2] + (-72.)*X.Yu*Yd2*Yu2Tr*z3*g2[2] + (72.)*X.Yu*Yu2*Yu2Tr*z3*g2[2]
                              + (-48.)*X.Yu*Yu2Yd2Tr*z3*g2[2] + (-72.)*X.Yu*Yu4Tr*z3*g2[2] + (77./60.)*X.Yu*Yd2*g2[0]*g2[2] + (907./60.)*X.Yu*Yu2*g2[0]*g2[2] + (-991./60.)*X.Yu*Yd2Tr*g2[0]*g2[2]
                              + (-2419./60.)*X.Yu*Yu2Tr*g2[0]*g2[2] + (-88./5.)*X.Yu*Yd2*z3*g2[0]*g2[2] + (-24./5.)*X.Yu*Yu2*z3*g2[0]*g2[2] + (12.)*X.Yu*Yd2Tr*z3*g2[0]*g2[2] + (204./5.)*X.Yu*Yu2Tr*z3*g2[0]*g2[2]
                              + (435./4.)*X.Yu*Yd2*g2[1]*g2[2] + (-183./4.)*X.Yu*Yu2*g2[1]*g2[2] + (-489./4.)*X.Yu*Yd2Tr*g2[1]*g2[2] + (-489./4.)*X.Yu*Yu2Tr*g2[1]*g2[2] + (-216.)*X.Yu*Yd2*z3*g2[1]*g2[2]
                              + (72.)*X.Yu*Yu2*z3*g2[1]*g2[2] + (108.)*X.Yu*Yd2Tr*z3*g2[1]*g2[2] + (108.)*X.Yu*Yu2Tr*z3*g2[1]*g2[2] + (-321./20.)*X.Yu*g2[0]*g2[1]*g2[2] + (-1089./800.)*X.La[1]*X.Yu*g4[0]
                              + (-363./400.)*X.La[2]*X.Yu*g4[0] + (-363./800.)*X.La[3]*X.Yu*g4[0] + (33223./19200.)*X.Yu*Yd2*g4[0] + (-251791./19200.)*X.Yu*Yu2*g4[0] + (-12423./3200.)*X.Yu*Yd2Tr*g4[0]
                              + (-1991./800.)*X.Yu*Ye2Tr*g4[0] + (-7731./640.)*X.Yu*Yu2Tr*g4[0] + (203./200.)*X.Yu*Yd2*z3*g4[0] + (-99./200.)*X.Yu*Yu2*z3*g4[0] + (-201./100.)*X.Yu*Yd2Tr*z3*g4[0]
                              + (3./100.)*X.Yu*Yu2Tr*z3*g4[0] + (2277./400.)*X.Yu*g2[1]*g4[0] + (-153./20.)*X.Yu*z3*g2[1]*g4[0] + (1597./120.)*X.Yu*g2[2]*g4[0] + (-748./25.)*X.Yu*z3*g2[2]*g4[0]
                              + (-171./32.)*X.La[1]*X.Yu*g4[1] + (-57./16.)*X.La[2]*X.Yu*g4[1] + (-57./32.)*X.La[3]*X.Yu*g4[1] + (3849./256.)*X.Yu*Yd2*g4[1] + (13287./256.)*X.Yu*Yu2*g4[1]
                              + (2997./128.)*X.Yu*Yd2Tr*g4[1] + (-105./32.)*X.Yu*Ye2Tr*g4[1] + (8757./128.)*X.Yu*Yu2Tr*g4[1] + (261./8.)*X.Yu*Yd2*z3*g4[1] + (-135./8.)*X.Yu*Yu2*z3*g4[1] + (-243./4.)*X.Yu*Yd2Tr*z3*g4[1]
                              + (-297./4.)*X.Yu*Yu2Tr*z3*g4[1] + (633./160.)*X.Yu*g2[0]*g4[1] + (-27./4.)*X.Yu*z3*g2[0]*g4[1] + (843./8.)*X.Yu*g2[2]*g4[1] + (-108.)*X.Yu*z3*g2[2]*g4[1] + (-2167./6.)*X.Yu*Yd2*g4[2]
                              + (2383./6.)*X.Yu*Yu2*g4[2] + (626./3.)*X.Yu*Yd2Tr*g4[2] + (722./3.)*X.Yu*Yu2Tr*g4[2] + (172.)*X.Yu*Yd2*z3*g4[2] + (-204.)*X.Yu*Yu2*z3*g4[2] + (-216.)*X.Yu*Yd2Tr*z3*g4[2]
                              + (-24.)*X.Yu*Yu2Tr*z3*g4[2] + (1633./60.)*X.Yu*g2[0]*g4[2] + (-176./5.)*X.Yu*z3*g2[0]*g4[2] + (987./4.)*X.Yu*g2[1]*g4[2] + (-144.)*X.Yu*z3*g2[1]*g4[2] + (598511./18000.)*X.Yu*g6[0]
                              + (-6613./500.)*X.Yu*z3*g6[0] + (307./32.)*X.Yu*g6[1] + (585./4.)*X.Yu*z3*g6[1] + (-4166./3.)*X.Yu*g6[2] + (640.)*X.Yu*z3*g6[2] + (-21./16.)*X.Yu*Yd2*La2[1] + (285./16.)*X.Yu*Yu2*La2[1]
                              + (-135./8.)*X.Yu*Yd2Tr*La2[1] + (-135./8.)*X.Yu*Yu2Tr*La2[1] + (9./4.)*X.Yu*g2[0]*La2[1] + (45./4.)*X.Yu*g2[1]*La2[1] + (-3./2.)*X.La[0]*X.Yu*La2[2] + (-3.)*X.La[1]*X.Yu*La2[2]
                              + (-3./2.)*X.La[3]*X.Yu*La2[2] + (-7./8.)*X.Yu*Yd2*La2[2] + (15./8.)*X.Yu*Yu2*La2[2] + (-15./4.)*X.Yu*Yd2Tr*La2[2] + (-5./2.)*X.Yu*Ye2Tr*La2[2] + (-15./4.)*X.Yu*Yu2Tr*La2[2]
                              + (3./2.)*X.Yu*g2[0]*La2[2] + (15./2.)*X.Yu*g2[1]*La2[2] + (-3./4.)*X.La[0]*X.Yu*La2[3] + (-3./2.)*X.La[1]*X.Yu*La2[3] + (-3.)*X.La[2]*X.Yu*La2[3] + (41./8.)*X.Yu*Yd2*La2[3]
                              + (-9./8.)*X.Yu*Yu2*La2[3] + (-15./4.)*X.Yu*Yd2Tr*La2[3] + (-5./2.)*X.Yu*Ye2Tr*La2[3] + (-15./4.)*X.Yu*Yu2Tr*La2[3] + (3./2.)*X.Yu*g2[0]*La2[3] + (15./2.)*X.Yu*g2[1]*La2[3]
                              + (-9./2.)*X.Yu*La3[1] + (-1.)*X.Yu*La3[2] + (-5./4.)*X.Yu*La3[3]);
        dX.Yd += loopfactor3*((-3./2.)*X.La[0]*X.La[2]*X.La[3]*X.Yd + (-3.)*X.La[1]*X.La[2]*X.La[3]*X.Yd + (-3./4.)*std::conj(X.La[4])*X.La[0]*X.La[4]*X.Yd
                              + (-3./2.)*std::conj(X.La[4])*X.La[1]*X.La[4]*X.Yd + (-9./2.)*std::conj(X.La[4])*X.La[2]*X.La[4]*X.Yd + (-27./4.)*std::conj(X.La[4])*X.La[3]*X.La[4]*X.Yd
                              + (15./8.)*X.La[2]*X.La[3]*X.Yd*Yd2 + (45./16.)*std::conj(X.La[4])*X.La[4]*X.Yd*Yd2 + (43./8.)*X.Yd*Yd2Yu2Yd2 + (-183./16.)*X.Yd*Yd2Yu4 + (63./2.)*X.La[1]*X.Yd*Yd4
                              + (83./16.)*X.Yd*Yd4Yu2 + (-345./16.)*X.Yd*Yd6 + (-7./8.)*X.La[2]*X.La[3]*X.Yd*Yu2 + (75./16.)*std::conj(X.La[4])*X.La[4]*X.Yd*Yu2 + (3./2.)*X.La[1]*X.Yd*Yu2Yd2
                              + (-37./8.)*X.Yd*Yu2Yd2Yu2 + (75./16.)*X.Yd*Yu2Yd4 + (-15.)*X.La[1]*X.Yd*Yu4 + (-95./8.)*X.Yd*Yu4Yd2 + (9./4.)*X.Yd*Yu6 + (-15./4.)*X.La[2]*X.La[3]*X.Yd*Yd2Tr
                              + (-45./8.)*std::conj(X.La[4])*X.La[4]*X.Yd*Yd2Tr + (45.)*X.La[1]*X.Yd*Yd2*Yd2Tr + (3./4.)*X.Yd*Yd2Yu2*Yd2Tr + (-9./4.)*X.Yd*Yd4*Yd2Tr + (21./4.)*X.Yd*Yu2Yd2*Yd2Tr
                              + (-69./4.)*X.Yd*Yu4*Yd2Tr + (-81./8.)*X.Yd*Yd2*Yd2Tr2 + (117./8.)*X.Yd*Yu2*Yd2Tr2 + (45./2.)*X.La[1]*X.Yd*Yd4Tr + (27.)*X.Yd*Yd2*Yd4Tr + (-153./8.)*X.Yd*Yu2*Yd4Tr
                              + (54.)*X.Yd*Yd2Tr*Yd4Tr + (-75./16.)*X.Yd*Yd6Tr + (-5./2.)*X.La[2]*X.La[3]*X.Yd*Ye2Tr + (-15./4.)*std::conj(X.La[4])*X.La[4]*X.Yd*Ye2Tr + (-15./4.)*X.La[2]*X.La[3]*X.Yd*Yu2Tr
                              + (-45./8.)*std::conj(X.La[4])*X.La[4]*X.Yd*Yu2Tr + (45.)*X.La[1]*X.Yd*Yd2*Yu2Tr + (3./4.)*X.Yd*Yd2Yu2*Yu2Tr + (-9./4.)*X.Yd*Yd4*Yu2Tr + (21./4.)*X.Yd*Yu2Yd2*Yu2Tr
                              + (-69./4.)*X.Yd*Yu4*Yu2Tr + (-81./4.)*X.Yd*Yd2*Yd2Tr*Yu2Tr + (117./4.)*X.Yd*Yu2*Yd2Tr*Yu2Tr + (54.)*X.Yd*Yd4Tr*Yu2Tr + (-81./8.)*X.Yd*Yd2*Yu2Tr2 + (117./8.)*X.Yd*Yu2*Yu2Tr2
                              + (177./4.)*X.Yd*Yu2*Yu2Yd2Tr + (-9./2.)*X.Yd*Yd2Tr*Yu2Yd2Tr + (-9./2.)*X.Yd*Yu2Tr*Yu2Yd2Tr + (39./16.)*X.Yd*Yu2Yd4Tr + (45./2.)*X.La[1]*X.Yd*Yu4Tr + (27.)*X.Yd*Yd2*Yu4Tr
                              + (-153./8.)*X.Yd*Yu2*Yu4Tr + (54.)*X.Yd*Yd2Tr*Yu4Tr + (54.)*X.Yd*Yu2Tr*Yu4Tr + (39./16.)*X.Yd*Yu4Yd2Tr + (-75./16.)*X.Yd*Yu6Tr + (12.)*X.Yd*Yd2Yu4*z3 + (9./2.)*X.Yd*Yd6*z3
                              + (12.)*X.Yd*Yu4Yd2*z3 + (-9./2.)*X.Yd*Yu6*z3 + (9.)*X.Yd*Yd6Tr*z3 + (-72.)*X.Yd*Yu2*Yu2Yd2Tr*z3 + (9.)*X.Yd*Yu6Tr*z3 + (3./2.)*X.La[2]*X.La[3]*X.Yd*g2[0]
                              + (9./4.)*std::conj(X.La[4])*X.La[4]*X.Yd*g2[0] + (-139./20.)*X.La[1]*X.Yd*Yd2*g2[0] + (63./40.)*X.Yd*Yd2Yu2*g2[0] + (-7./20.)*X.Yd*Yd4*g2[0] + (137./160.)*X.Yd*Yu2Yd2*g2[0]
                              + (1043./160.)*X.Yd*Yu4*g2[0] + (-9.)*X.Yd*Yd2*Yd2Tr*g2[0] + (-83./40.)*X.Yd*Yu2*Yd2Tr*g2[0] + (-81./20.)*X.Yd*Yd2Tr2*g2[0] + (-909./80.)*X.Yd*Yd4Tr*g2[0] + (-81./5.)*X.Yd*Yd2*Yu2Tr*g2[0]
                              + (127./40.)*X.Yd*Yu2*Yu2Tr*g2[0] + (-81./10.)*X.Yd*Yd2Tr*Yu2Tr*g2[0] + (-81./20.)*X.Yd*Yu2Tr2*g2[0] + (-93./40.)*X.Yd*Yu2Yd2Tr*g2[0] + (-633./80.)*X.Yd*Yu4Tr*g2[0]
                              + (-31./20.)*X.Yd*Yd2Yu2*z3*g2[0] + (-11./4.)*X.Yd*Yu2Yd2*z3*g2[0] + (-17./10.)*X.Yd*Yu4*z3*g2[0] + (-27./5.)*X.Yd*Yd2*Yd2Tr*z3*g2[0] + (36./5.)*X.Yd*Yu2*Yd2Tr*z3*g2[0]
                              + (27./5.)*X.Yd*Yd4Tr*z3*g2[0] + (9./5.)*X.Yd*Yu2*Yu2Tr*z3*g2[0] + (24./5.)*X.Yd*Yu2Yd2Tr*z3*g2[0] + (-9./5.)*X.Yd*Yu4Tr*z3*g2[0] + (15./2.)*X.La[2]*X.La[3]*X.Yd*g2[1]
                              + (45./4.)*std::conj(X.La[4])*X.La[4]*X.Yd*g2[1] + (-135./4.)*X.La[1]*X.Yd*Yd2*g2[1] + (195./16.)*X.Yd*Yd2Yu2*g2[1] + (-27./4.)*X.Yd*Yd4*g2[1] + (39./32.)*X.Yd*Yu2Yd2*g2[1]
                              + (579./32.)*X.Yd*Yu4*g2[1] + (-81./4.)*X.Yd*Yd2*Yd2Tr*g2[1] + (-135./8.)*X.Yd*Yu2*Yd2Tr*g2[1] + (-81./4.)*X.Yd*Yd2Tr2*g2[1] + (-837./16.)*X.Yd*Yd4Tr*g2[1] + (-81./4.)*X.Yd*Yd2*Yu2Tr*g2[1]
                              + (-135./8.)*X.Yd*Yu2*Yu2Tr*g2[1] + (-81./2.)*X.Yd*Yd2Tr*Yu2Tr*g2[1] + (-81./4.)*X.Yd*Yu2Tr2*g2[1] + (-63./8.)*X.Yd*Yu2Yd2Tr*g2[1] + (-837./16.)*X.Yd*Yu4Tr*g2[1]
                              + (-9./4.)*X.Yd*Yd2Yu2*z3*g2[1] + (-9./4.)*X.Yd*Yu2Yd2*z3*g2[1] + (-45./2.)*X.Yd*Yu4*z3*g2[1] + (-27.)*X.Yd*Yd2*Yd2Tr*z3*g2[1] + (27.)*X.Yd*Yu2*Yd2Tr*z3*g2[1] + (27.)*X.Yd*Yd4Tr*z3*g2[1]
                              + (-27.)*X.Yd*Yd2*Yu2Tr*z3*g2[1] + (27.)*X.Yd*Yu2*Yu2Tr*z3*g2[1] + (27.)*X.Yd*Yu4Tr*z3*g2[1] + (-27./80.)*X.La[1]*X.Yd*g2[0]*g2[1] + (-27./80.)*X.La[3]*X.Yd*g2[0]*g2[1]
                              + (-141./640.)*X.Yd*Yd2*g2[0]*g2[1] + (669./128.)*X.Yd*Yu2*g2[0]*g2[1] + (4317./320.)*X.Yd*Yd2Tr*g2[0]*g2[1] + (-39./320.)*X.Yd*Yu2Tr*g2[0]*g2[1] + (18./5.)*X.Yd*Yd2*z3*g2[0]*g2[1]
                              + (-81./20.)*X.Yd*Yu2*z3*g2[0]*g2[1] + (-27./5.)*X.Yd*Yd2Tr*z3*g2[0]*g2[1] + (27./2.)*X.Yd*Yu2Tr*z3*g2[0]*g2[1] + (8.)*X.La[1]*X.Yd*Yd2*g2[2] + (-18.)*X.Yd*Yd2Yu2*g2[2]
                              + (-76.)*X.Yd*Yd4*g2[2] + (28.)*X.Yd*Yu2Yd2*g2[2] + (26.)*X.Yd*Yu4*g2[2] + (-177./2.)*X.Yd*Yd2*Yd2Tr*g2[2] + (97./2.)*X.Yd*Yu2*Yd2Tr*g2[2] + (15./2.)*X.Yd*Yd4Tr*g2[2]
                              + (-177./2.)*X.Yd*Yd2*Yu2Tr*g2[2] + (97./2.)*X.Yd*Yu2*Yu2Tr*g2[2] + (57.)*X.Yd*Yu2Yd2Tr*g2[2] + (15./2.)*X.Yd*Yu4Tr*g2[2] + (8.)*X.Yd*Yd2Yu2*z3*g2[2] + (8.)*X.Yd*Yu2Yd2*z3*g2[2]
                              + (80.)*X.Yd*Yu4*z3*g2[2] + (72.)*X.Yd*Yd2*Yd2Tr*z3*g2[2] + (-72.)*X.Yd*Yu2*Yd2Tr*z3*g2[2] + (-72.)*X.Yd*Yd4Tr*z3*g2[2] + (72.)*X.Yd*Yd2*Yu2Tr*z3*g2[2] + (-72.)*X.Yd*Yu2*Yu2Tr*z3*g2[2]
                              + (-48.)*X.Yd*Yu2Yd2Tr*z3*g2[2] + (-72.)*X.Yd*Yu4Tr*z3*g2[2] + (-89./60.)*X.Yd*Yd2*g2[0]*g2[2] + (809./60.)*X.Yd*Yu2*g2[0]*g2[2] + (-991./60.)*X.Yd*Yd2Tr*g2[0]*g2[2]
                              + (-2419./60.)*X.Yd*Yu2Tr*g2[0]*g2[2] + (72./5.)*X.Yd*Yd2*z3*g2[0]*g2[2] + (-184./5.)*X.Yd*Yu2*z3*g2[0]*g2[2] + (12.)*X.Yd*Yd2Tr*z3*g2[0]*g2[2] + (204./5.)*X.Yd*Yu2Tr*z3*g2[0]*g2[2]
                              + (-183./4.)*X.Yd*Yd2*g2[1]*g2[2] + (435./4.)*X.Yd*Yu2*g2[1]*g2[2] + (-489./4.)*X.Yd*Yd2Tr*g2[1]*g2[2] + (-489./4.)*X.Yd*Yu2Tr*g2[1]*g2[2] + (72.)*X.Yd*Yd2*z3*g2[1]*g2[2]
                              + (-216.)*X.Yd*Yu2*z3*g2[1]*g2[2] + (108.)*X.Yd*Yd2Tr*z3*g2[1]*g2[2] + (108.)*X.Yd*Yu2Tr*z3*g2[1]*g2[2] + (-153./20.)*X.Yd*g2[0]*g2[1]*g2[2] + (-9./32.)*X.La[1]*X.Yd*g4[0]
                              + (-3./16.)*X.La[2]*X.Yd*g4[0] + (-3./32.)*X.La[3]*X.Yd*g4[0] + (-38123./3840.)*X.Yd*Yd2*g4[0] + (5815./768.)*X.Yd*Yu2*g4[0] + (-207./128.)*X.Yd*Yd2Tr*g4[0]
                              + (-79./160.)*X.Yd*Ye2Tr*g4[0] + (-42607./3200.)*X.Yd*Yu2Tr*g4[0] + (3./200.)*X.Yd*Yd2*z3*g4[0] + (209./200.)*X.Yd*Yu2*z3*g4[0] + (-87./100.)*X.Yd*Yd2Tr*z3*g4[0]
                              + (-111./100.)*X.Yd*Yu2Tr*z3*g4[0] + (2931./400.)*X.Yd*g2[1]*g4[0] + (-9./4.)*X.Yd*z3*g2[1]*g4[0] + (-2959./600.)*X.Yd*g2[2]*g4[0] + (-44./5.)*X.Yd*z3*g2[2]*g4[0]
                              + (-171./32.)*X.La[1]*X.Yd*g4[1] + (-57./16.)*X.La[2]*X.Yd*g4[1] + (-57./32.)*X.La[3]*X.Yd*g4[1] + (13287./256.)*X.Yd*Yd2*g4[1] + (3849./256.)*X.Yd*Yu2*g4[1]
                              + (8757./128.)*X.Yd*Yd2Tr*g4[1] + (-105./32.)*X.Yd*Ye2Tr*g4[1] + (2997./128.)*X.Yd*Yu2Tr*g4[1] + (-135./8.)*X.Yd*Yd2*z3*g4[1] + (261./8.)*X.Yd*Yu2*z3*g4[1] + (-297./4.)*X.Yd*Yd2Tr*z3*g4[1]
                              + (-243./4.)*X.Yd*Yu2Tr*z3*g4[1] + (-9./160.)*X.Yd*g2[0]*g4[1] + (-27./4.)*X.Yd*z3*g2[0]*g4[1] + (843./8.)*X.Yd*g2[2]*g4[1] + (-108.)*X.Yd*z3*g2[2]*g4[1] + (2383./6.)*X.Yd*Yd2*g4[2]
                              + (-2167./6.)*X.Yd*Yu2*g4[2] + (722./3.)*X.Yd*Yd2Tr*g4[2] + (626./3.)*X.Yd*Yu2Tr*g4[2] + (-204.)*X.Yd*Yd2*z3*g4[2] + (172.)*X.Yd*Yu2*z3*g4[2] + (-24.)*X.Yd*Yd2Tr*z3*g4[2]
                              + (-216.)*X.Yd*Yu2Tr*z3*g4[2] + (833./12.)*X.Yd*g2[0]*g4[2] + (-176./5.)*X.Yd*z3*g2[0]*g4[2] + (987./4.)*X.Yd*g2[1]*g4[2] + (-144.)*X.Yd*z3*g2[1]*g4[2] + (110431./9000.)*X.Yd*g6[0]
                              + (-389./100.)*X.Yd*z3*g6[0] + (307./32.)*X.Yd*g6[1] + (585./4.)*X.Yd*z3*g6[1] + (-4166./3.)*X.Yd*g6[2] + (640.)*X.Yd*z3*g6[2] + (285./16.)*X.Yd*Yd2*La2[1] + (-21./16.)*X.Yd*Yu2*La2[1]
                              + (-135./8.)*X.Yd*Yd2Tr*La2[1] + (-135./8.)*X.Yd*Yu2Tr*La2[1] + (9./4.)*X.Yd*g2[0]*La2[1] + (45./4.)*X.Yd*g2[1]*La2[1] + (-3./2.)*X.La[0]*X.Yd*La2[2] + (-3.)*X.La[1]*X.Yd*La2[2]
                              + (-3./2.)*X.La[3]*X.Yd*La2[2] + (15./8.)*X.Yd*Yd2*La2[2] + (-7./8.)*X.Yd*Yu2*La2[2] + (-15./4.)*X.Yd*Yd2Tr*La2[2] + (-5./2.)*X.Yd*Ye2Tr*La2[2] + (-15./4.)*X.Yd*Yu2Tr*La2[2]
                              + (3./2.)*X.Yd*g2[0]*La2[2] + (15./2.)*X.Yd*g2[1]*La2[2] + (-3./4.)*X.La[0]*X.Yd*La2[3] + (-3./2.)*X.La[1]*X.Yd*La2[3] + (-3.)*X.La[2]*X.Yd*La2[3] + (-9./8.)*X.Yd*Yd2*La2[3]
                              + (41./8.)*X.Yd*Yu2*La2[3] + (-15./4.)*X.Yd*Yd2Tr*La2[3] + (-5./2.)*X.Yd*Ye2Tr*La2[3] + (-15./4.)*X.Yd*Yu2Tr*La2[3] + (3./2.)*X.Yd*g2[0]*La2[3] + (15./2.)*X.Yd*g2[1]*La2[3]
                              + (-9./2.)*X.Yd*La3[1] + (-1.)*X.Yd*La3[2] + (-5./4.)*X.Yd*La3[3]);
        dX.Ye += loopfactor3*((-3.)*X.La[0]*X.La[2]*X.La[3]*X.Ye + (-3./2.)*X.La[1]*X.La[2]*X.La[3]*X.Ye + (-3./2.)*std::conj(X.La[4])*X.La[0]*X.La[4]*X.Ye
                              + (-3./4.)*std::conj(X.La[4])*X.La[1]*X.La[4]*X.Ye + (-9./2.)*std::conj(X.La[4])*X.La[2]*X.La[4]*X.Ye + (-27./4.)*std::conj(X.La[4])*X.La[3]*X.La[4]*X.Ye
                              + (15./8.)*X.La[2]*X.La[3]*X.Ye*Ye2 + (45./16.)*std::conj(X.La[4])*X.La[4]*X.Ye*Ye2 + (63./2.)*X.La[0]*X.Ye*Ye4 + (-345./16.)*X.Ye*Ye6 + (-15./2.)*X.La[2]*X.La[3]*X.Ye*Yd2Tr
                              + (-45./4.)*std::conj(X.La[4])*X.La[4]*X.Ye*Yd2Tr + (-5./4.)*X.La[2]*X.La[3]*X.Ye*Ye2Tr + (-15./8.)*std::conj(X.La[4])*X.La[4]*X.Ye*Ye2Tr + (15.)*X.La[0]*X.Ye*Ye2*Ye2Tr
                              + (-3./4.)*X.Ye*Ye4*Ye2Tr + (-9./8.)*X.Ye*Ye2*Ye2Tr2 + (15./2.)*X.La[0]*X.Ye*Ye4Tr + (9.)*X.Ye*Ye2*Ye4Tr + (6.)*X.Ye*Ye2Tr*Ye4Tr + (-25./16.)*X.Ye*Ye6Tr
                              + (-15./2.)*X.La[2]*X.La[3]*X.Ye*Yu2Tr + (-45./4.)*std::conj(X.La[4])*X.La[4]*X.Ye*Yu2Tr + (9./2.)*X.Ye*Ye6*z3 + (3.)*X.Ye*Ye6Tr*z3 + (3./2.)*X.La[2]*X.La[3]*X.Ye*g2[0]
                              + (9./4.)*std::conj(X.La[4])*X.La[4]*X.Ye*g2[0] + (-99./20.)*X.La[0]*X.Ye*Ye2*g2[0] + (-369./20.)*X.Ye*Ye4*g2[0] + (-171./20.)*X.Ye*Ye2*Ye2Tr*g2[0] + (-9./20.)*X.Ye*Ye2Tr2*g2[0]
                              + (-99./80.)*X.Ye*Ye4Tr*g2[0] + (27./5.)*X.Ye*Ye2*Ye2Tr*z3*g2[0] + (-27./5.)*X.Ye*Ye4Tr*z3*g2[0] + (15./2.)*X.La[2]*X.La[3]*X.Ye*g2[1] + (45./4.)*std::conj(X.La[4])*X.La[4]*X.Ye*g2[1]
                              + (-135./4.)*X.La[0]*X.Ye*Ye2*g2[1] + (-27./4.)*X.Ye*Ye4*g2[1] + (-27./4.)*X.Ye*Ye2*Ye2Tr*g2[1] + (-9./4.)*X.Ye*Ye2Tr2*g2[1] + (-279./16.)*X.Ye*Ye4Tr*g2[1] + (-9.)*X.Ye*Ye2*Ye2Tr*z3*g2[1]
                              + (9.)*X.Ye*Ye4Tr*z3*g2[1] + (261./80.)*X.La[0]*X.Ye*g2[0]*g2[1] + (261./80.)*X.La[3]*X.Ye*g2[0]*g2[1] + (-7173./640.)*X.Ye*Ye2*g2[0]*g2[1] + (2223./320.)*X.Ye*Ye2Tr*g2[0]*g2[1]
                              + (243./10.)*X.Ye*Ye2*z3*g2[0]*g2[1] + (54./5.)*X.Ye*Ye2Tr*z3*g2[0]*g2[1] + (-621./160.)*X.La[0]*X.Ye*g4[0] + (-207./80.)*X.La[2]*X.Ye*g4[0] + (-207./160.)*X.La[3]*X.Ye*g4[0]
                              + (-22761./1280.)*X.Ye*Ye2*g4[0] + (-813./160.)*X.Ye*Yd2Tr*g4[0] + (-6777./640.)*X.Ye*Ye2Tr*g4[0] + (-7989./800.)*X.Ye*Yu2Tr*g4[0] + (-1377./200.)*X.Ye*Ye2*z3*g4[0]
                              + (351./100.)*X.Ye*Ye2Tr*z3*g4[0] + (3699./400.)*X.Ye*g2[1]*g4[0] + (-81./4.)*X.Ye*z3*g2[1]*g4[0] + (7227./100.)*X.Ye*g2[2]*g4[0] + (-396./5.)*X.Ye*z3*g2[2]*g4[0]
                              + (-171./32.)*X.La[0]*X.Ye*g4[1] + (-57./16.)*X.La[2]*X.Ye*g4[1] + (-57./32.)*X.La[3]*X.Ye*g4[1] + (13287./256.)*X.Ye*Ye2*g4[1] + (-315./32.)*X.Ye*Yd2Tr*g4[1]
                              + (2919./128.)*X.Ye*Ye2Tr*g4[1] + (-315./32.)*X.Ye*Yu2Tr*g4[1] + (-135./8.)*X.Ye*Ye2*z3*g4[1] + (-99./4.)*X.Ye*Ye2Tr*z3*g4[1] + (1557./160.)*X.Ye*g2[0]*g4[1] + (-27./4.)*X.Ye*z3*g2[0]*g4[1]
                              + (351./4.)*X.Ye*g2[2]*g4[1] + (-108.)*X.Ye*z3*g2[2]*g4[1] + (158643./2000.)*X.Ye*g6[0] + (-3501./100.)*X.Ye*z3*g6[0] + (307./32.)*X.Ye*g6[1] + (585./4.)*X.Ye*z3*g6[1]
                              + (285./16.)*X.Ye*Ye2*La2[0] + (-45./8.)*X.Ye*Ye2Tr*La2[0] + (9./4.)*X.Ye*g2[0]*La2[0] + (45./4.)*X.Ye*g2[1]*La2[0] + (-3.)*X.La[0]*X.Ye*La2[2] + (-3./2.)*X.La[1]*X.Ye*La2[2]
                              + (-3./2.)*X.La[3]*X.Ye*La2[2] + (15./8.)*X.Ye*Ye2*La2[2] + (-15./2.)*X.Ye*Yd2Tr*La2[2] + (-5./4.)*X.Ye*Ye2Tr*La2[2] + (-15./2.)*X.Ye*Yu2Tr*La2[2] + (3./2.)*X.Ye*g2[0]*La2[2]
                              + (15./2.)*X.Ye*g2[1]*La2[2] + (-3./2.)*X.La[0]*X.Ye*La2[3] + (-3./4.)*X.La[1]*X.Ye*La2[3] + (-3.)*X.La[2]*X.Ye*La2[3] + (-9./8.)*X.Ye*Ye2*La2[3] + (-15./2.)*X.Ye*Yd2Tr*La2[3]
                              + (-5./4.)*X.Ye*Ye2Tr*La2[3] + (-15./2.)*X.Ye*Yu2Tr*La2[3] + (3./2.)*X.Ye*g2[0]*La2[3] + (15./2.)*X.Ye*g2[1]*La2[3] + (-9./2.)*X.Ye*La3[0] + (-1.)*X.Ye*La3[2] + (-5./4.)*X.Ye*La3[3]);
        }
      }
    }
  } else {
    dX.setZero();
  }
};
