#include "sm.h"
#include "ckm.h"
#include "pmns.h"

#define loopfactor 0.006332573977646111 // 1/(4 Pi)^2
#define loopfactor2 0.00004010149318236069 // 1/(4 Pi)^4
#define loopfactor3 0.00000025394567219137 // 1/(4 Pi)^6
#define loopfactor4 0.00000000160812975545 // 1/(4 Pi)^8
#define z3 1.2020569031595942854 // Zeta(3)

using namespace Eigen;

// beta functions for sm
void sm::operator()(const sm &X, sm &dX, const double) {
  if(check()) {

    // variables for powers of parameters for nloops >= 1
    gauge<3> g2(X.g.array().square().matrix());
    gauge<3> g3(X.g.array().cube().matrix());
    gauge<3> g4(g2.array().square().matrix());
    std::complex<double> La2 = X.La[0]*X.La[0];
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
    
    if(!X.weylordering || (X.weylordering && X.nloops > 1)) {
    dX.Yu = loopfactor*((-3./2.)*X.Yu*Yd2 + (3./2.)*X.Yu*Yu2 + (3.)*X.Yu*Yd2Tr + X.Yu*Ye2Tr + (3.)*X.Yu*Yu2Tr + (-17./20.)*X.Yu*g2[0] + (-9./4.)*X.Yu*g2[1] + (-8.)*X.Yu*g2[2]);
    dX.Yd = loopfactor*((3./2.)*X.Yd*Yd2 + (-3./2.)*X.Yd*Yu2 + (3.)*X.Yd*Yd2Tr + X.Yd*Ye2Tr + (3.)*X.Yd*Yu2Tr + (-1./4.)*X.Yd*g2[0] + (-9./4.)*X.Yd*g2[1] + (-8.)*X.Yd*g2[2]);
    dX.Ye = loopfactor*((3./2.)*X.Ye*Ye2 + (3.)*X.Ye*Yd2Tr + X.Ye*Ye2Tr + (3.)*X.Ye*Yu2Tr + (-9./4.)*X.Ye*g2[0] + (-9./4.)*X.Ye*g2[1]);
    }
    
    if(!X.weylordering || (X.weylordering && X.nloops > 2)) {
    dX.La[0] = loopfactor*((12.)*X.La[0]*Yd2Tr + (-24.)*Yd4Tr + (4.)*X.La[0]*Ye2Tr + (-8.)*Ye4Tr + (12.)*X.La[0]*Yu2Tr + (-24.)*Yu4Tr + (-9./5.)*X.La[0]*g2[0] + (-9.)*X.La[0]*g2[1]
                         + (9./5.)*g2[0]*g2[1] + (27./50.)*g4[0] + (9./2.)*g4[1] + (6.)*La2);
    }

    if(X.nloops > 1) {
      gauge<3> g5(g3.cwiseProduct(g2));
      gauge<3> g6(g3.array().square().matrix());
      std::complex<double> La3 = X.La[0]*La2;
      yukawa Yu6 = Yu4*Yu2;             std::complex<double> Yu6Tr = Yu6.trace();
      yukawa Yd6 = Yd4*Yd2;             std::complex<double> Yd6Tr = Yd6.trace();
      yukawa Ye6 = Ye4*Ye2;             std::complex<double> Ye6Tr = Ye6.trace();
      yukawa Yu4Yd2 = Yu4*Yd2;          std::complex<double> Yu4Yd2Tr = Yu4Yd2.trace();
      yukawa Yu2Yd4 = Yu2*Yd4;          std::complex<double> Yu2Yd4Tr = Yu2Yd4.trace();

      dX.g[0] += loopfactor2*((-1./2.)*Yd2Tr*g3[0] + (-3./2.)*Ye2Tr*g3[0] + (-17./10.)*Yu2Tr*g3[0] + (27./10.)*g2[1]*g3[0] + (44./5.)*g2[2]*g3[0] + (199./50.)*g5[0]).real();
      dX.g[1] += loopfactor2*((-3./2.)*Yd2Tr*g3[1] + (-1./2.)*Ye2Tr*g3[1] + (-3./2.)*Yu2Tr*g3[1] + (9./10.)*g2[0]*g3[1] + (12.)*g2[2]*g3[1] + (35./6.)*g5[1]).real();
      dX.g[2] += loopfactor2*((-2.)*Yd2Tr*g3[2] + (-2.)*Yu2Tr*g3[2] + (11./10.)*g2[0]*g3[2] + (9./2.)*g2[1]*g3[2] + (-26.)*g5[2]).real();
      
      if(!X.weylordering || (X.weylordering && X.nloops > 2)) {
      dX.Yu += loopfactor2*((-1./4.)*X.Yu*Yd2Yu2 + (11./4.)*X.Yu*Yd4 + (-6.)*X.La[0]*X.Yu*Yu2 + (-1.)*X.Yu*Yu2Yd2 + (3./2.)*X.Yu*Yu4 + (15./4.)*X.Yu*Yd2*Yd2Tr + (-27./4.)*X.Yu*Yu2*Yd2Tr
                          + (-27./4.)*X.Yu*Yd4Tr + (5./4.)*X.Yu*Yd2*Ye2Tr + (-9./4.)*X.Yu*Yu2*Ye2Tr + (-9./4.)*X.Yu*Ye4Tr + (15./4.)*X.Yu*Yd2*Yu2Tr + (-27./4.)*X.Yu*Yu2*Yu2Tr + (3./2.)*X.Yu*Yu2Yd2Tr
                          + (-27./4.)*X.Yu*Yu4Tr + (-43./80.)*X.Yu*Yd2*g2[0] + (223./80.)*X.Yu*Yu2*g2[0] + (5./8.)*X.Yu*Yd2Tr*g2[0] + (15./8.)*X.Yu*Ye2Tr*g2[0] + (17./8.)*X.Yu*Yu2Tr*g2[0]
                          + (9./16.)*X.Yu*Yd2*g2[1] + (135./16.)*X.Yu*Yu2*g2[1] + (45./8.)*X.Yu*Yd2Tr*g2[1] + (15./8.)*X.Yu*Ye2Tr*g2[1] + (45./8.)*X.Yu*Yu2Tr*g2[1] + (-9./20.)*X.Yu*g2[0]*g2[1]
                          + (-16.)*X.Yu*Yd2*g2[2] + (16.)*X.Yu*Yu2*g2[2] + (20.)*X.Yu*Yd2Tr*g2[2] + (20.)*X.Yu*Yu2Tr*g2[2] + (19./15.)*X.Yu*g2[0]*g2[2] + (9.)*X.Yu*g2[1]*g2[2] + (1187./600.)*X.Yu*g4[0]
                          + (-23./4.)*X.Yu*g4[1] + (-108.)*X.Yu*g4[2] + (3./2.)*X.Yu*La2);
      dX.Yd += loopfactor2*((-6.)*X.La[0]*X.Yd*Yd2 + (-1.)*X.Yd*Yd2Yu2 + (3./2.)*X.Yd*Yd4 + (-1./4.)*X.Yd*Yu2Yd2 + (11./4.)*X.Yd*Yu4 + (-27./4.)*X.Yd*Yd2*Yd2Tr + (15./4.)*X.Yd*Yu2*Yd2Tr
                          + (-27./4.)*X.Yd*Yd4Tr + (-9./4.)*X.Yd*Yd2*Ye2Tr + (5./4.)*X.Yd*Yu2*Ye2Tr + (-9./4.)*X.Yd*Ye4Tr + (-27./4.)*X.Yd*Yd2*Yu2Tr + (15./4.)*X.Yd*Yu2*Yu2Tr + (3./2.)*X.Yd*Yu2Yd2Tr
                          + (-27./4.)*X.Yd*Yu4Tr + (187./80.)*X.Yd*Yd2*g2[0] + (-79./80.)*X.Yd*Yu2*g2[0] + (5./8.)*X.Yd*Yd2Tr*g2[0] + (15./8.)*X.Yd*Ye2Tr*g2[0] + (17./8.)*X.Yd*Yu2Tr*g2[0]
                          + (135./16.)*X.Yd*Yd2*g2[1] + (9./16.)*X.Yd*Yu2*g2[1] + (45./8.)*X.Yd*Yd2Tr*g2[1] + (15./8.)*X.Yd*Ye2Tr*g2[1] + (45./8.)*X.Yd*Yu2Tr*g2[1] + (-27./20.)*X.Yd*g2[0]*g2[1]
                          + (16.)*X.Yd*Yd2*g2[2] + (-16.)*X.Yd*Yu2*g2[2] + (20.)*X.Yd*Yd2Tr*g2[2] + (20.)*X.Yd*Yu2Tr*g2[2] + (31./15.)*X.Yd*g2[0]*g2[2] + (9.)*X.Yd*g2[1]*g2[2] + (-127./600.)*X.Yd*g4[0]
                          + (-23./4.)*X.Yd*g4[1] + (-108.)*X.Yd*g4[2] + (3./2.)*X.Yd*La2);
      dX.Ye += loopfactor2*((-6.)*X.La[0]*X.Ye*Ye2 + (3./2.)*X.Ye*Ye4 + (-27./4.)*X.Ye*Ye2*Yd2Tr + (-27./4.)*X.Ye*Yd4Tr + (-9./4.)*X.Ye*Ye2*Ye2Tr + (-9./4.)*X.Ye*Ye4Tr + (-27./4.)*X.Ye*Ye2*Yu2Tr
                          + (3./2.)*X.Ye*Yu2Yd2Tr + (-27./4.)*X.Ye*Yu4Tr + (387./80.)*X.Ye*Ye2*g2[0] + (5./8.)*X.Ye*Yd2Tr*g2[0] + (15./8.)*X.Ye*Ye2Tr*g2[0] + (17./8.)*X.Ye*Yu2Tr*g2[0]
                          + (135./16.)*X.Ye*Ye2*g2[1] + (45./8.)*X.Ye*Yd2Tr*g2[1] + (15./8.)*X.Ye*Ye2Tr*g2[1] + (45./8.)*X.Ye*Yu2Tr*g2[1] + (27./20.)*X.Ye*g2[0]*g2[1] + (20.)*X.Ye*Yd2Tr*g2[2]
                          + (20.)*X.Ye*Yu2Tr*g2[2] + (1371./200.)*X.Ye*g4[0] + (-23./4.)*X.Ye*g4[1] + (3./2.)*X.Ye*La2);
      }
      
      if(!X.weylordering || (X.weylordering && X.nloops > 3)) {
      dX.La[0] += loopfactor2*((-3.)*X.La[0]*Yd4Tr + (120.)*Yd6Tr + (-1.)*X.La[0]*Ye4Tr + (40.)*Ye6Tr + (-42.)*X.La[0]*Yu2Yd2Tr + (-24.)*Yu2Yd4Tr + (-3.)*X.La[0]*Yu4Tr + (-24.)*Yu4Yd2Tr + (120.)*Yu6Tr
                             + (5./2.)*X.La[0]*Yd2Tr*g2[0] + (16./5.)*Yd4Tr*g2[0] + (15./2.)*X.La[0]*Ye2Tr*g2[0] + (-48./5.)*Ye4Tr*g2[0] + (17./2.)*X.La[0]*Yu2Tr*g2[0] + (-32./5.)*Yu4Tr*g2[0]
                             + (45./2.)*X.La[0]*Yd2Tr*g2[1] + (15./2.)*X.La[0]*Ye2Tr*g2[1] + (45./2.)*X.La[0]*Yu2Tr*g2[1] + (117./20.)*X.La[0]*g2[0]*g2[1] + (54./5.)*Yd2Tr*g2[0]*g2[1] + (66./5.)*Ye2Tr*g2[0]*g2[1]
                             + (126./5.)*Yu2Tr*g2[0]*g2[1] + (80.)*X.La[0]*Yd2Tr*g2[2] + (-128.)*Yd4Tr*g2[2] + (80.)*X.La[0]*Yu2Tr*g2[2] + (-128.)*Yu4Tr*g2[2] + (1887./200.)*X.La[0]*g4[0] + (9./5.)*Yd2Tr*g4[0]
                             + (-9.)*Ye2Tr*g4[0] + (-171./25.)*Yu2Tr*g4[0] + (-1677./100.)*g2[1]*g4[0] + (-73./8.)*X.La[0]*g4[1] + (-9.)*Yd2Tr*g4[1] + (-3.)*Ye2Tr*g4[1] + (-9.)*Yu2Tr*g4[1] + (-289./20.)*g2[0]*g4[1]
                             + (-3411./500.)*g6[0] + (305./4.)*g6[1] + (-36.)*Yd2Tr*La2 + (-12.)*Ye2Tr*La2 + (-36.)*Yu2Tr*La2 + (27./5.)*g2[0]*La2 + (27.)*g2[1]*La2 + (-39./2.)*La3);
      }

      if(X.nloops > 2) {
        gauge<3> g7(g5.cwiseProduct(g2));
        gauge<3> g8(g4.array().square().matrix());
        std::complex<double> La4 = La2*La2;
        std::complex<double> Yu2Tr2 = Yu2Tr*Yu2Tr;
        std::complex<double> Yd2Tr2 = Yd2Tr*Yd2Tr;
        std::complex<double> Ye2Tr2 = Ye2Tr*Ye2Tr;
        std::complex<double> Yu4Tr2 = Yu4Tr*Yu4Tr;
        std::complex<double> Yd4Tr2 = Yd4Tr*Yd4Tr;
        std::complex<double> Ye4Tr2 = Ye4Tr*Ye4Tr;
        std::complex<double> Yu2Yd2Tr2 = Yu2Yd2Tr*Yu2Yd2Tr;
        yukawa Yd2Yu4 = Yd2*Yu4;                yukawa Yu2Yd2Yu2 = Yu2*Yd2*Yu2;
        yukawa Yd4Yu2 = Yd4*Yu2;                yukawa Yd2Yu2Yd2 = Yd2*Yu2*Yd2;
        yukawa Yu8 = Yu6*Yu2;                   std::complex<double> Yu8Tr = Yu8.trace();
        yukawa Yd8 = Yd6*Yd2;                   std::complex<double> Yd8Tr = Yd8.trace();
        yukawa Ye8 = Ye6*Ye2;                   std::complex<double> Ye8Tr = Ye8.trace();
        yukawa Yu4Yd4 = Yu4*Yd4;                std::complex<double> Yu4Yd4Tr = Yu4Yd4.trace();
        yukawa Yu2Yd6 = Yu2*Yd6;                std::complex<double> Yu2Yd6Tr = Yu2Yd6.trace();
        yukawa Yu6Yd2 = Yu6*Yd2;                std::complex<double> Yu6Yd2Tr = Yu6Yd2.trace();
        yukawa Yu2Yd2Yu2Yd2 = Yu2*Yd2*Yu2*Yd2;  std::complex<double> Yu2Yd2Yu2Yd2Tr = Yu2Yd2Yu2Yd2.trace();

        dX.g[0] += loopfactor3*((51./40.)*Yd2Tr2*g3[0] + (183./80.)*Yd4Tr*g3[0] + (157./20.)*Yd2Tr*Ye2Tr*g3[0] + (99./40.)*Ye2Tr2*g3[0] + (261./80.)*Ye4Tr*g3[0] + (177./20.)*Yd2Tr*Yu2Tr*g3[0]
                              + (199./20.)*Ye2Tr*Yu2Tr*g3[0] + (303./40.)*Yu2Tr2*g3[0] + (3./8.)*Yu2Yd2Tr*g3[0] + (339./80.)*Yu4Tr*g3[0] + (9./20.)*X.La[0]*g2[1]*g3[0] + (-1311./160.)*Yd2Tr*g2[1]*g3[0]
                              + (-1629./160.)*Ye2Tr*g2[1]*g3[0] + (-471./32.)*Yu2Tr*g2[1]*g3[0] + (-17./5.)*Yd2Tr*g2[2]*g3[0] + (-29./5.)*Yu2Tr*g2[2]*g3[0] + (-3./5.)*g2[1]*g2[2]*g3[0]
                              + (789./64.)*g3[0]*g4[1] + (297./5.)*g3[0]*g4[2] + (27./100.)*X.La[0]*g5[0] +  (-1267./800.)*Yd2Tr*g5[0] + (-2529./800.)*Ye2Tr*g5[0] + (-2827./800.)*Yu2Tr*g5[0]
                              + (123./160.)*g2[1]*g5[0] + (-137./75.)*g2[2]*g5[0] + (-388613./24000.)*g7[0] + (-9./20.)*g3[0]*La2).real();
        dX.g[1] += loopfactor3*((45./8.)*Yd2Tr2*g3[1] + (57./16.)*Yd4Tr*g3[1] + (15./4.)*Yd2Tr*Ye2Tr*g3[1] + (5./8.)*Ye2Tr2*g3[1] + (19./16.)*Ye4Tr*g3[1] + (45./4.)*Yd2Tr*Yu2Tr*g3[1]
                              + (15./4.)*Ye2Tr*Yu2Tr*g3[1] + (45./8.)*Yu2Tr2*g3[1] + (27./8.)*Yu2Yd2Tr*g3[1] + (57./16.)*Yu4Tr*g3[1] + (3./20.)*X.La[0]*g2[0]*g3[1] + (-533./160.)*Yd2Tr*g2[0]*g3[1]
                              + (-51./32.)*Ye2Tr*g2[0]*g3[1] + (-593./160.)*Yu2Tr*g2[0]*g3[1] + (-7.)*Yd2Tr*g2[2]*g3[1] + (-7.)*Yu2Tr*g2[2]*g3[1] + (-1./5.)*g2[0]*g2[2]*g3[1] + (-5597./1600.)*g3[1]*g4[0]
                              + (81.)*g3[1]*g4[2] + (3./4.)*X.La[0]*g5[1] +  (-729./32.)*Yd2Tr*g5[1] + (-243./32.)*Ye2Tr*g5[1] + (-729./32.)*Yu2Tr*g5[1] + (873./160.)*g2[0]*g5[1] + (39.)*g2[2]*g5[1]
                              + (324953./1728.)*g7[1] + (-3./4.)*g3[1]*La2).real();
        dX.g[2] += loopfactor3*((21./2.)*Yd2Tr2*g3[2] + (9./2.)*Yd4Tr*g3[2] + (7./2.)*Yd2Tr*Ye2Tr*g3[2] + (21.)*Yd2Tr*Yu2Tr*g3[2] + (7./2.)*Ye2Tr*Yu2Tr*g3[2] + (21./2.)*Yu2Tr2*g3[2] + (-3.)*Yu2Yd2Tr*g3[2]
                              + (9./2.)*Yu4Tr*g3[2] + (-89./40.)*Yd2Tr*g2[0]*g3[2] + (-101./40.)*Yu2Tr*g2[0]*g3[2] + (-93./8.)*Yd2Tr*g2[1]*g3[2] + (-93./8.)*Yu2Tr*g2[1]*g3[2] + (-3./40.)*g2[0]*g2[1]*g3[2]
                              + (-523./120.)*g3[2]*g4[0] + (109./8.)*g3[2]*g4[1] + (-40.)*Yd2Tr*g5[2] + (-40.)*Yu2Tr*g5[2] + (77./15.)*g2[0]*g5[2] + (21.)*g2[1]*g5[2] + (65./2.)*g7[2]).real();
                              
        if(!X.weylordering || (X.weylordering && X.nloops > 3)) {
        dX.Yu += loopfactor3*((3./2.)*X.La[0]*X.Yu*Yd2Yu2 + (-37./8.)*X.Yu*Yd2Yu2Yd2 + (75./16.)*X.Yu*Yd2Yu4 + (-15.)*X.La[0]*X.Yu*Yd4 + (-95./8.)*X.Yu*Yd4Yu2 + (9./4.)*X.Yu*Yd6 + (43./8.)*X.Yu*Yu2Yd2Yu2
                            + (-183./16.)*X.Yu*Yu2Yd4 + (63./2.)*X.La[0]*X.Yu*Yu4 + (83./16.)*X.Yu*Yu4Yd2 + (-345./16.)*X.Yu*Yu6 + (21./4.)*X.Yu*Yd2Yu2*Yd2Tr + (-69./4.)*X.Yu*Yd4*Yd2Tr
                            + (45.)*X.La[0]*X.Yu*Yu2*Yd2Tr + (3./4.)*X.Yu*Yu2Yd2*Yd2Tr + (-9./4.)*X.Yu*Yu4*Yd2Tr + (117./8.)*X.Yu*Yd2*Yd2Tr2 + (-81./8.)*X.Yu*Yu2*Yd2Tr2 + (45./2.)*X.La[0]*X.Yu*Yd4Tr
                            + (-153./8.)*X.Yu*Yd2*Yd4Tr + (27.)*X.Yu*Yu2*Yd4Tr + (54.)*X.Yu*Yd2Tr*Yd4Tr + (-75./16.)*X.Yu*Yd6Tr + (7./4.)*X.Yu*Yd2Yu2*Ye2Tr + (-23./4.)*X.Yu*Yd4*Ye2Tr
                            + (15.)*X.La[0]*X.Yu*Yu2*Ye2Tr + (1./4.)*X.Yu*Yu2Yd2*Ye2Tr + (-3./4.)*X.Yu*Yu4*Ye2Tr + (39./4.)*X.Yu*Yd2*Yd2Tr*Ye2Tr + (-27./4.)*X.Yu*Yu2*Yd2Tr*Ye2Tr + (18.)*X.Yu*Yd4Tr*Ye2Tr
                            + (13./8.)*X.Yu*Yd2*Ye2Tr2 + (-9./8.)*X.Yu*Yu2*Ye2Tr2 + (15./2.)*X.La[0]*X.Yu*Ye4Tr + (-51./8.)*X.Yu*Yd2*Ye4Tr + (9.)*X.Yu*Yu2*Ye4Tr + (18.)*X.Yu*Yd2Tr*Ye4Tr + (6.)*X.Yu*Ye2Tr*Ye4Tr
                            + (-25./16.)*X.Yu*Ye6Tr + (39./16.)*X.Yu*Yu2Yd4Tr + (21./4.)*X.Yu*Yd2Yu2*Yu2Tr + (-69./4.)*X.Yu*Yd4*Yu2Tr + (45.)*X.La[0]*X.Yu*Yu2*Yu2Tr + (3./4.)*X.Yu*Yu2Yd2*Yu2Tr
                            + (-9./4.)*X.Yu*Yu4*Yu2Tr + (117./4.)*X.Yu*Yd2*Yd2Tr*Yu2Tr + (-81./4.)*X.Yu*Yu2*Yd2Tr*Yu2Tr + (54.)*X.Yu*Yd4Tr*Yu2Tr + (39./4.)*X.Yu*Yd2*Ye2Tr*Yu2Tr + (-27./4.)*X.Yu*Yu2*Ye2Tr*Yu2Tr
                            + (18.)*X.Yu*Ye4Tr*Yu2Tr + (117./8.)*X.Yu*Yd2*Yu2Tr2 + (-81./8.)*X.Yu*Yu2*Yu2Tr2 + (177./4.)*X.Yu*Yd2*Yu2Yd2Tr + (-9./2.)*X.Yu*Yd2Tr*Yu2Yd2Tr + (-3./2.)*X.Yu*Ye2Tr*Yu2Yd2Tr
                            + (-9./2.)*X.Yu*Yu2Tr*Yu2Yd2Tr + (45./2.)*X.La[0]*X.Yu*Yu4Tr + (-153./8.)*X.Yu*Yd2*Yu4Tr + (27.)*X.Yu*Yu2*Yu4Tr + (54.)*X.Yu*Yd2Tr*Yu4Tr + (18.)*X.Yu*Ye2Tr*Yu4Tr
                            + (54.)*X.Yu*Yu2Tr*Yu4Tr + (39./16.)*X.Yu*Yu4Yd2Tr + (-75./16.)*X.Yu*Yu6Tr + (12.)*X.Yu*Yd4Yu2*z3 + (-9./2.)*X.Yu*Yd6*z3 + (12.)*X.Yu*Yu2Yd4*z3 + (9./2.)*X.Yu*Yu6*z3
                            + (9.)*X.Yu*Yd6Tr*z3 + (3.)*X.Yu*Ye6Tr*z3 + (-72.)*X.Yu*Yd2*Yu2Yd2Tr*z3 + (9.)*X.Yu*Yu6Tr*z3 + (-199./160.)*X.Yu*Yd2Yu2*g2[0] + (1047./160.)*X.Yu*Yd4*g2[0]
                            + (-127./20.)*X.La[0]*X.Yu*Yu2*g2[0] + (3./5.)*X.Yu*Yu2Yd2*g2[0] + (-28./5.)*X.Yu*Yu4*g2[0] + (23./8.)*X.Yu*Yd2*Yd2Tr*g2[0] + (-57./10.)*X.Yu*Yu2*Yd2Tr*g2[0]
                            + (-81./20.)*X.Yu*Yd2Tr2*g2[0] + (-909./80.)*X.Yu*Yd4Tr*g2[0] + (163./24.)*X.Yu*Yd2*Ye2Tr*g2[0] + (-99./10.)*X.Yu*Yu2*Ye2Tr*g2[0] + (-27./10.)*X.Yu*Yd2Tr*Ye2Tr*g2[0]
                            + (-9./20.)*X.Yu*Ye2Tr2*g2[0] + (-99./80.)*X.Yu*Ye4Tr*g2[0] + (65./8.)*X.Yu*Yd2*Yu2Tr*g2[0] + (-129./10.)*X.Yu*Yu2*Yu2Tr*g2[0] + (-81./10.)*X.Yu*Yd2Tr*Yu2Tr*g2[0]
                            + (-27./10.)*X.Yu*Ye2Tr*Yu2Tr*g2[0] + (-81./20.)*X.Yu*Yu2Tr2*g2[0] + (-93./40.)*X.Yu*Yu2Yd2Tr*g2[0] + (-633./80.)*X.Yu*Yu4Tr*g2[0] + (41./20.)*X.Yu*Yd2Yu2*z3*g2[0]
                            + (-53./10.)*X.Yu*Yd4*z3*g2[0] + (17./20.)*X.Yu*Yu2Yd2*z3*g2[0] + (9./5.)*X.Yu*Yd2*Yd2Tr*z3*g2[0] + (-18./5.)*X.Yu*Yu2*Yd2Tr*z3*g2[0] + (27./5.)*X.Yu*Yd4Tr*z3*g2[0]
                            + (-27./5.)*X.Yu*Yd2*Ye2Tr*z3*g2[0] + (24./5.)*X.Yu*Yu2*Ye2Tr*z3*g2[0] + (-27./5.)*X.Yu*Ye4Tr*z3*g2[0] + (-18./5.)*X.Yu*Yd2*Yu2Tr*z3*g2[0] + (9./5.)*X.Yu*Yu2*Yu2Tr*z3*g2[0]
                            + (24./5.)*X.Yu*Yu2Yd2Tr*z3*g2[0] + (-9./5.)*X.Yu*Yu4Tr*z3*g2[0] + (39./32.)*X.Yu*Yd2Yu2*g2[1] + (579./32.)*X.Yu*Yd4*g2[1] + (-135./4.)*X.La[0]*X.Yu*Yu2*g2[1]
                            + (195./16.)*X.Yu*Yu2Yd2*g2[1] + (-27./4.)*X.Yu*Yu4*g2[1] + (-135./8.)*X.Yu*Yd2*Yd2Tr*g2[1] + (-81./4.)*X.Yu*Yu2*Yd2Tr*g2[1] + (-81./4.)*X.Yu*Yd2Tr2*g2[1] + (-837./16.)*X.Yu*Yd4Tr*g2[1]
                            + (-45./8.)*X.Yu*Yd2*Ye2Tr*g2[1] + (-27./4.)*X.Yu*Yu2*Ye2Tr*g2[1] + (-27./2.)*X.Yu*Yd2Tr*Ye2Tr*g2[1] + (-9./4.)*X.Yu*Ye2Tr2*g2[1] + (-279./16.)*X.Yu*Ye4Tr*g2[1]
                            + (-135./8.)*X.Yu*Yd2*Yu2Tr*g2[1] + (-81./4.)*X.Yu*Yu2*Yu2Tr*g2[1] + (-81./2.)*X.Yu*Yd2Tr*Yu2Tr*g2[1] + (-27./2.)*X.Yu*Ye2Tr*Yu2Tr*g2[1] + (-81./4.)*X.Yu*Yu2Tr2*g2[1]
                            + (-63./8.)*X.Yu*Yu2Yd2Tr*g2[1] + (-837./16.)*X.Yu*Yu4Tr*g2[1] + (-9./4.)*X.Yu*Yd2Yu2*z3*g2[1] + (-45./2.)*X.Yu*Yd4*z3*g2[1] + (-9./4.)*X.Yu*Yu2Yd2*z3*g2[1] + (27.)*X.Yu*Yd2*Yd2Tr*z3*g2[1]
                            + (-27.)*X.Yu*Yu2*Yd2Tr*z3*g2[1] + (27.)*X.Yu*Yd4Tr*z3*g2[1] + (9.)*X.Yu*Yd2*Ye2Tr*z3*g2[1] + (-9.)*X.Yu*Yu2*Ye2Tr*z3*g2[1] + (9.)*X.Yu*Ye4Tr*z3*g2[1] + (27.)*X.Yu*Yd2*Yu2Tr*z3*g2[1]
                            + (-27.)*X.Yu*Yu2*Yu2Tr*z3*g2[1] + (27.)*X.Yu*Yu4Tr*z3*g2[1] + (117./80.)*X.La[0]*X.Yu*g2[0]*g2[1] + (-1443./640.)*X.Yu*Yd2*g2[0]*g2[1] + (-2193./640.)*X.Yu*Yu2*g2[0]*g2[1]
                            + (2589./320.)*X.Yu*Yd2Tr*g2[0]*g2[1] + (-1041./320.)*X.Yu*Ye2Tr*g2[0]*g2[1] + (1029./64.)*X.Yu*Yu2Tr*g2[0]*g2[1] + (27./10.)*X.Yu*Yd2*z3*g2[0]*g2[1] + (207./20.)*X.Yu*Yu2*z3*g2[0]*g2[1]
                            + (-9./5.)*X.Yu*Ye2Tr*z3*g2[0]*g2[1] + (81./10.)*X.Yu*Yu2Tr*z3*g2[0]*g2[1] + (28.)*X.Yu*Yd2Yu2*g2[2] + (26.)*X.Yu*Yd4*g2[2] + (8.)*X.La[0]*X.Yu*Yu2*g2[2] + (-18.)*X.Yu*Yu2Yd2*g2[2]
                            + (-76.)*X.Yu*Yu4*g2[2] + (97./2.)*X.Yu*Yd2*Yd2Tr*g2[2] + (-177./2.)*X.Yu*Yu2*Yd2Tr*g2[2] + (15./2.)*X.Yu*Yd4Tr*g2[2] + (-43./6.)*X.Yu*Yd2*Ye2Tr*g2[2] + (5./2.)*X.Yu*Yu2*Ye2Tr*g2[2]
                            + (97./2.)*X.Yu*Yd2*Yu2Tr*g2[2] + (-177./2.)*X.Yu*Yu2*Yu2Tr*g2[2] + (57.)*X.Yu*Yu2Yd2Tr*g2[2] + (15./2.)*X.Yu*Yu4Tr*g2[2] + (8.)*X.Yu*Yd2Yu2*z3*g2[2] + (80.)*X.Yu*Yd4*z3*g2[2]
                            + (8.)*X.Yu*Yu2Yd2*z3*g2[2] + (-72.)*X.Yu*Yd2*Yd2Tr*z3*g2[2] + (72.)*X.Yu*Yu2*Yd2Tr*z3*g2[2] + (-72.)*X.Yu*Yd4Tr*z3*g2[2] + (-72.)*X.Yu*Yd2*Yu2Tr*z3*g2[2] + (72.)*X.Yu*Yu2*Yu2Tr*z3*g2[2]
                            + (-48.)*X.Yu*Yu2Yd2Tr*z3*g2[2] + (-72.)*X.Yu*Yu4Tr*z3*g2[2] + (77./60.)*X.Yu*Yd2*g2[0]*g2[2] + (907./60.)*X.Yu*Yu2*g2[0]*g2[2] + (-991./60.)*X.Yu*Yd2Tr*g2[0]*g2[2]
                            + (-2419./60.)*X.Yu*Yu2Tr*g2[0]*g2[2] + (-88./5.)*X.Yu*Yd2*z3*g2[0]*g2[2] + (-24./5.)*X.Yu*Yu2*z3*g2[0]*g2[2] + (12.)*X.Yu*Yd2Tr*z3*g2[0]*g2[2] + (204./5.)*X.Yu*Yu2Tr*z3*g2[0]*g2[2]
                            + (435./4.)*X.Yu*Yd2*g2[1]*g2[2] + (-183./4.)*X.Yu*Yu2*g2[1]*g2[2] + (-489./4.)*X.Yu*Yd2Tr*g2[1]*g2[2] + (-489./4.)*X.Yu*Yu2Tr*g2[1]*g2[2] + (-216.)*X.Yu*Yd2*z3*g2[1]*g2[2]
                            + (72.)*X.Yu*Yu2*z3*g2[1]*g2[2] + (108.)*X.Yu*Yd2Tr*z3*g2[1]*g2[2] + (108.)*X.Yu*Yu2Tr*z3*g2[1]*g2[2] + (-321./20.)*X.Yu*g2[0]*g2[1]*g2[2] + (-1089./800.)*X.La[0]*X.Yu*g4[0]
                            + (30877./19200.)*X.Yu*Yd2*g4[0] + (-238741./19200.)*X.Yu*Yu2*g4[0] + (-477./128.)*X.Yu*Yd2Tr*g4[0] + (-9659./640.)*X.Yu*Ye2Tr*g4[0] + (-36573./3200.)*X.Yu*Yu2Tr*g4[0]
                            + (203./200.)*X.Yu*Yd2*z3*g4[0] + (-99./200.)*X.Yu*Yu2*z3*g4[0] + (-201./100.)*X.Yu*Yd2Tr*z3*g4[0] + (-807./100.)*X.Yu*Ye2Tr*z3*g4[0] + (3./100.)*X.Yu*Yu2Tr*z3*g4[0]
                            + (1227./320.)*X.Yu*g2[1]*g4[0] + (-1377./200.)*X.Yu*z3*g2[1]*g4[0] + (2047./150.)*X.Yu*g2[2]*g4[0] + (-748./25.)*X.Yu*z3*g2[2]*g4[0] + (-171./32.)*X.La[0]*X.Yu*g4[1]
                            + (3663./256.)*X.Yu*Yd2*g4[1] + (14193./256.)*X.Yu*Yu2*g4[1] + (3339./128.)*X.Yu*Yd2Tr*g4[1] + (1113./128.)*X.Yu*Ye2Tr*g4[1] + (9099./128.)*X.Yu*Yu2Tr*g4[1] + (261./8.)*X.Yu*Yd2*z3*g4[1]
                            + (-135./8.)*X.Yu*Yu2*z3*g4[1] + (-243./4.)*X.Yu*Yd2Tr*z3*g4[1] + (-81./4.)*X.Yu*Ye2Tr*z3*g4[1] + (-297./4.)*X.Yu*Yu2Tr*z3*g4[1] + (819./320.)*X.Yu*g2[0]*g4[1]
                            + (-243./40.)*X.Yu*z3*g2[0]*g4[1] + (435./4.)*X.Yu*g2[2]*g4[1] + (-108.)*X.Yu*z3*g2[2]*g4[1] + (-2167./6.)*X.Yu*Yd2*g4[2] + (2383./6.)*X.Yu*Yu2*g4[2] + (626./3.)*X.Yu*Yd2Tr*g4[2]
                            + (722./3.)*X.Yu*Yu2Tr*g4[2] + (172.)*X.Yu*Yd2*z3*g4[2] + (-204.)*X.Yu*Yu2*z3*g4[2] + (-216.)*X.Yu*Yd2Tr*z3*g4[2] + (-24.)*X.Yu*Yu2Tr*z3*g4[2] + (1633./60.)*X.Yu*g2[0]*g4[2]
                            + (-176./5.)*X.Yu*z3*g2[0]*g4[2] + (987./4.)*X.Yu*g2[1]*g4[2] + (-144.)*X.Yu*z3*g2[1]*g4[2] + (763523./24000.)*X.Yu*g6[0] + (-13073./1000.)*X.Yu*z3*g6[0] + (455./576.)*X.Yu*g6[1]
                            + (1125./8.)*X.Yu*z3*g6[1] + (-4166./3.)*X.Yu*g6[2] + (640.)*X.Yu*z3*g6[2] + (-21./16.)*X.Yu*Yd2*La2 + (285./16.)*X.Yu*Yu2*La2 + (-135./8.)*X.Yu*Yd2Tr*La2
                            + (-45./8.)*X.Yu*Ye2Tr*La2 + (-135./8.)*X.Yu*Yu2Tr*La2 + (9./4.)*X.Yu*g2[0]*La2 + (45./4.)*X.Yu*g2[1]*La2 + (-9./2.)*X.Yu*La3);
        dX.Yd += loopfactor3*((43./8.)*X.Yd*Yd2Yu2Yd2 + (-183./16.)*X.Yd*Yd2Yu4 + (63./2.)*X.La[0]*X.Yd*Yd4 + (83./16.)*X.Yd*Yd4Yu2 + (-345./16.)*X.Yd*Yd6 + (3./2.)*X.La[0]*X.Yd*Yu2Yd2
                            + (-37./8.)*X.Yd*Yu2Yd2Yu2 + (75./16.)*X.Yd*Yu2Yd4 + (-15.)*X.La[0]*X.Yd*Yu4 + (-95./8.)*X.Yd*Yu4Yd2 + (9./4.)*X.Yd*Yu6 + (45.)*X.La[0]*X.Yd*Yd2*Yd2Tr + (3./4.)*X.Yd*Yd2Yu2*Yd2Tr
                            + (-9./4.)*X.Yd*Yd4*Yd2Tr + (21./4.)*X.Yd*Yu2Yd2*Yd2Tr + (-69./4.)*X.Yd*Yu4*Yd2Tr + (-81./8.)*X.Yd*Yd2*Yd2Tr2 + (117./8.)*X.Yd*Yu2*Yd2Tr2 + (45./2.)*X.La[0]*X.Yd*Yd4Tr
                            + (27.)*X.Yd*Yd2*Yd4Tr + (-153./8.)*X.Yd*Yu2*Yd4Tr + (54.)*X.Yd*Yd2Tr*Yd4Tr + (-75./16.)*X.Yd*Yd6Tr + (15.)*X.La[0]*X.Yd*Yd2*Ye2Tr + (1./4.)*X.Yd*Yd2Yu2*Ye2Tr
                            + (-3./4.)*X.Yd*Yd4*Ye2Tr + (7./4.)*X.Yd*Yu2Yd2*Ye2Tr + (-23./4.)*X.Yd*Yu4*Ye2Tr + (-27./4.)*X.Yd*Yd2*Yd2Tr*Ye2Tr + (39./4.)*X.Yd*Yu2*Yd2Tr*Ye2Tr + (18.)*X.Yd*Yd4Tr*Ye2Tr
                            + (-9./8.)*X.Yd*Yd2*Ye2Tr2 + (13./8.)*X.Yd*Yu2*Ye2Tr2 + (15./2.)*X.La[0]*X.Yd*Ye4Tr + (9.)*X.Yd*Yd2*Ye4Tr + (-51./8.)*X.Yd*Yu2*Ye4Tr + (18.)*X.Yd*Yd2Tr*Ye4Tr + (6.)*X.Yd*Ye2Tr*Ye4Tr
                            + (-25./16.)*X.Yd*Ye6Tr + (39./16.)*X.Yd*Yu2Yd4Tr + (45.)*X.La[0]*X.Yd*Yd2*Yu2Tr + (3./4.)*X.Yd*Yd2Yu2*Yu2Tr + (-9./4.)*X.Yd*Yd4*Yu2Tr + (21./4.)*X.Yd*Yu2Yd2*Yu2Tr
                            + (-69./4.)*X.Yd*Yu4*Yu2Tr + (-81./4.)*X.Yd*Yd2*Yd2Tr*Yu2Tr + (117./4.)*X.Yd*Yu2*Yd2Tr*Yu2Tr + (54.)*X.Yd*Yd4Tr*Yu2Tr + (-27./4.)*X.Yd*Yd2*Ye2Tr*Yu2Tr + (39./4.)*X.Yd*Yu2*Ye2Tr*Yu2Tr
                            + (18.)*X.Yd*Ye4Tr*Yu2Tr + (-81./8.)*X.Yd*Yd2*Yu2Tr2 + (117./8.)*X.Yd*Yu2*Yu2Tr2 + (177./4.)*X.Yd*Yu2*Yu2Yd2Tr + (-9./2.)*X.Yd*Yd2Tr*Yu2Yd2Tr + (-3./2.)*X.Yd*Ye2Tr*Yu2Yd2Tr
                            + (-9./2.)*X.Yd*Yu2Tr*Yu2Yd2Tr + (45./2.)*X.La[0]*X.Yd*Yu4Tr + (27.)*X.Yd*Yd2*Yu4Tr + (-153./8.)*X.Yd*Yu2*Yu4Tr + (54.)*X.Yd*Yd2Tr*Yu4Tr + (18.)*X.Yd*Ye2Tr*Yu4Tr
                            + (54.)*X.Yd*Yu2Tr*Yu4Tr + (39./16.)*X.Yd*Yu4Yd2Tr + (-75./16.)*X.Yd*Yu6Tr + (12.)*X.Yd*Yd2Yu4*z3 + (9./2.)*X.Yd*Yd6*z3 + (12.)*X.Yd*Yu4Yd2*z3 + (-9./2.)*X.Yd*Yu6*z3
                            + (9.)*X.Yd*Yd6Tr*z3 + (3.)*X.Yd*Ye6Tr*z3 + (-72.)*X.Yd*Yu2*Yu2Yd2Tr*z3 + (9.)*X.Yd*Yu6Tr*z3 + (-139./20.)*X.La[0]*X.Yd*Yd2*g2[0] + (63./40.)*X.Yd*Yd2Yu2*g2[0]
                            + (-7./20.)*X.Yd*Yd4*g2[0] + (137./160.)*X.Yd*Yu2Yd2*g2[0] + (1043./160.)*X.Yd*Yu4*g2[0] + (-9.)*X.Yd*Yd2*Yd2Tr*g2[0] + (-83./40.)*X.Yd*Yu2*Yd2Tr*g2[0] + (-81./20.)*X.Yd*Yd2Tr2*g2[0]
                            + (-909./80.)*X.Yd*Yd4Tr*g2[0] + (-11.)*X.Yd*Yd2*Ye2Tr*g2[0] + (617./120.)*X.Yd*Yu2*Ye2Tr*g2[0] + (-27./10.)*X.Yd*Yd2Tr*Ye2Tr*g2[0] + (-9./20.)*X.Yd*Ye2Tr2*g2[0]
                            + (-99./80.)*X.Yd*Ye4Tr*g2[0] + (-81./5.)*X.Yd*Yd2*Yu2Tr*g2[0] + (127./40.)*X.Yd*Yu2*Yu2Tr*g2[0] + (-81./10.)*X.Yd*Yd2Tr*Yu2Tr*g2[0] + (-27./10.)*X.Yd*Ye2Tr*Yu2Tr*g2[0]
                            + (-81./20.)*X.Yd*Yu2Tr2*g2[0] + (-93./40.)*X.Yd*Yu2Yd2Tr*g2[0] + (-633./80.)*X.Yd*Yu4Tr*g2[0] + (-31./20.)*X.Yd*Yd2Yu2*z3*g2[0] + (-11./4.)*X.Yd*Yu2Yd2*z3*g2[0]
                            + (-17./10.)*X.Yd*Yu4*z3*g2[0] + (-27./5.)*X.Yd*Yd2*Yd2Tr*z3*g2[0] + (36./5.)*X.Yd*Yu2*Yd2Tr*z3*g2[0] + (27./5.)*X.Yd*Yd4Tr*z3*g2[0] + (21./5.)*X.Yd*Yd2*Ye2Tr*z3*g2[0]
                            + (-18./5.)*X.Yd*Yu2*Ye2Tr*z3*g2[0] + (-27./5.)*X.Yd*Ye4Tr*z3*g2[0] + (9./5.)*X.Yd*Yu2*Yu2Tr*z3*g2[0] + (24./5.)*X.Yd*Yu2Yd2Tr*z3*g2[0] + (-9./5.)*X.Yd*Yu4Tr*z3*g2[0]
                            + (-135./4.)*X.La[0]*X.Yd*Yd2*g2[1] + (195./16.)*X.Yd*Yd2Yu2*g2[1] + (-27./4.)*X.Yd*Yd4*g2[1] + (39./32.)*X.Yd*Yu2Yd2*g2[1] + (579./32.)*X.Yd*Yu4*g2[1] + (-81./4.)*X.Yd*Yd2*Yd2Tr*g2[1]
                            + (-135./8.)*X.Yd*Yu2*Yd2Tr*g2[1] + (-81./4.)*X.Yd*Yd2Tr2*g2[1] + (-837./16.)*X.Yd*Yd4Tr*g2[1] + (-27./4.)*X.Yd*Yd2*Ye2Tr*g2[1] + (-45./8.)*X.Yd*Yu2*Ye2Tr*g2[1]
                            + (-27./2.)*X.Yd*Yd2Tr*Ye2Tr*g2[1] + (-9./4.)*X.Yd*Ye2Tr2*g2[1] + (-279./16.)*X.Yd*Ye4Tr*g2[1] + (-81./4.)*X.Yd*Yd2*Yu2Tr*g2[1] + (-135./8.)*X.Yd*Yu2*Yu2Tr*g2[1]
                            + (-81./2.)*X.Yd*Yd2Tr*Yu2Tr*g2[1] + (-27./2.)*X.Yd*Ye2Tr*Yu2Tr*g2[1] + (-81./4.)*X.Yd*Yu2Tr2*g2[1] + (-63./8.)*X.Yd*Yu2Yd2Tr*g2[1] + (-837./16.)*X.Yd*Yu4Tr*g2[1]
                            + (-9./4.)*X.Yd*Yd2Yu2*z3*g2[1] + (-9./4.)*X.Yd*Yu2Yd2*z3*g2[1] + (-45./2.)*X.Yd*Yu4*z3*g2[1] + (-27.)*X.Yd*Yd2*Yd2Tr*z3*g2[1] + (27.)*X.Yd*Yu2*Yd2Tr*z3*g2[1] + (27.)*X.Yd*Yd4Tr*z3*g2[1]
                            + (-9.)*X.Yd*Yd2*Ye2Tr*z3*g2[1] + (9.)*X.Yd*Yu2*Ye2Tr*z3*g2[1] + (9.)*X.Yd*Ye4Tr*z3*g2[1] + (-27.)*X.Yd*Yd2*Yu2Tr*z3*g2[1] + (27.)*X.Yd*Yu2*Yu2Tr*z3*g2[1] + (27.)*X.Yd*Yu4Tr*z3*g2[1]
                            + (-27./80.)*X.La[0]*X.Yd*g2[0]*g2[1] + (-141./640.)*X.Yd*Yd2*g2[0]*g2[1] + (669./128.)*X.Yd*Yu2*g2[0]*g2[1] + (4317./320.)*X.Yd*Yd2Tr*g2[0]*g2[1] + (-1233./320.)*X.Yd*Ye2Tr*g2[0]*g2[1]
                            + (-39./320.)*X.Yd*Yu2Tr*g2[0]*g2[1] + (18./5.)*X.Yd*Yd2*z3*g2[0]*g2[1] + (-81./20.)*X.Yd*Yu2*z3*g2[0]*g2[1] + (-27./5.)*X.Yd*Yd2Tr*z3*g2[0]*g2[1] + (54./5.)*X.Yd*Ye2Tr*z3*g2[0]*g2[1]
                            + (27./2.)*X.Yd*Yu2Tr*z3*g2[0]*g2[1] + (8.)*X.La[0]*X.Yd*Yd2*g2[2] + (-18.)*X.Yd*Yd2Yu2*g2[2] + (-76.)*X.Yd*Yd4*g2[2] + (28.)*X.Yd*Yu2Yd2*g2[2] + (26.)*X.Yd*Yu4*g2[2]
                            + (-177./2.)*X.Yd*Yd2*Yd2Tr*g2[2] + (97./2.)*X.Yd*Yu2*Yd2Tr*g2[2] + (15./2.)*X.Yd*Yd4Tr*g2[2] + (5./2.)*X.Yd*Yd2*Ye2Tr*g2[2] + (-43./6.)*X.Yd*Yu2*Ye2Tr*g2[2]
                            + (-177./2.)*X.Yd*Yd2*Yu2Tr*g2[2] + (97./2.)*X.Yd*Yu2*Yu2Tr*g2[2] + (57.)*X.Yd*Yu2Yd2Tr*g2[2] + (15./2.)*X.Yd*Yu4Tr*g2[2] + (8.)*X.Yd*Yd2Yu2*z3*g2[2] + (8.)*X.Yd*Yu2Yd2*z3*g2[2]
                            + (80.)*X.Yd*Yu4*z3*g2[2] + (72.)*X.Yd*Yd2*Yd2Tr*z3*g2[2] + (-72.)*X.Yd*Yu2*Yd2Tr*z3*g2[2] + (-72.)*X.Yd*Yd4Tr*z3*g2[2] + (72.)*X.Yd*Yd2*Yu2Tr*z3*g2[2] + (-72.)*X.Yd*Yu2*Yu2Tr*z3*g2[2]
                            + (-48.)*X.Yd*Yu2Yd2Tr*z3*g2[2] + (-72.)*X.Yd*Yu4Tr*z3*g2[2] + (-89./60.)*X.Yd*Yd2*g2[0]*g2[2] + (809./60.)*X.Yd*Yu2*g2[0]*g2[2] + (-991./60.)*X.Yd*Yd2Tr*g2[0]*g2[2]
                            + (-2419./60.)*X.Yd*Yu2Tr*g2[0]*g2[2] + (72./5.)*X.Yd*Yd2*z3*g2[0]*g2[2] + (-184./5.)*X.Yd*Yu2*z3*g2[0]*g2[2] + (12.)*X.Yd*Yd2Tr*z3*g2[0]*g2[2] + (204./5.)*X.Yd*Yu2Tr*z3*g2[0]*g2[2]
                            + (-183./4.)*X.Yd*Yd2*g2[1]*g2[2] + (435./4.)*X.Yd*Yu2*g2[1]*g2[2] + (-489./4.)*X.Yd*Yd2Tr*g2[1]*g2[2] + (-489./4.)*X.Yd*Yu2Tr*g2[1]*g2[2] + (72.)*X.Yd*Yd2*z3*g2[1]*g2[2]
                            + (-216.)*X.Yd*Yu2*z3*g2[1]*g2[2] + (108.)*X.Yd*Yd2Tr*z3*g2[1]*g2[2] + (108.)*X.Yd*Yu2Tr*z3*g2[1]*g2[2] + (-153./20.)*X.Yd*g2[0]*g2[1]*g2[2] + (-9./32.)*X.La[0]*X.Yd*g4[0]
                            + (-181597./19200.)*X.Yd*Yd2*g4[0] + (138421./19200.)*X.Yd*Yu2*g4[0] + (-4677./3200.)*X.Yd*Yd2Tr*g4[0] + (-45239./3200.)*X.Yd*Ye2Tr*g4[0] + (-1621./128.)*X.Yd*Yu2Tr*g4[0]
                            + (3./200.)*X.Yd*Yd2*z3*g4[0] + (209./200.)*X.Yd*Yu2*z3*g4[0] + (-87./100.)*X.Yd*Yd2Tr*z3*g4[0] + (351./100.)*X.Yd*Ye2Tr*z3*g4[0] + (-111./100.)*X.Yd*Yu2Tr*z3*g4[0]
                            + (2139./320.)*X.Yd*g2[1]*g4[0] + (-81./40.)*X.Yd*z3*g2[1]*g4[0] + (-337./75.)*X.Yd*g2[2]*g4[0] + (-44./5.)*X.Yd*z3*g2[2]*g4[0] + (-171./32.)*X.La[0]*X.Yd*g4[1]
                            + (14193./256.)*X.Yd*Yd2*g4[1] + (3663./256.)*X.Yd*Yu2*g4[1] + (9099./128.)*X.Yd*Yd2Tr*g4[1] + (3033./128.)*X.Yd*Ye2Tr*g4[1] + (3339./128.)*X.Yd*Yu2Tr*g4[1] + (-135./8.)*X.Yd*Yd2*z3*g4[1]
                            + (261./8.)*X.Yd*Yu2*z3*g4[1] + (-297./4.)*X.Yd*Yd2Tr*z3*g4[1] + (-99./4.)*X.Yd*Ye2Tr*z3*g4[1] + (-243./4.)*X.Yd*Yu2Tr*z3*g4[1] + (-633./320.)*X.Yd*g2[0]*g4[1]
                            + (-243./40.)*X.Yd*z3*g2[0]*g4[1] + (435./4.)*X.Yd*g2[2]*g4[1] + (-108.)*X.Yd*z3*g2[2]*g4[1] + (2383./6.)*X.Yd*Yd2*g4[2] + (-2167./6.)*X.Yd*Yu2*g4[2] + (722./3.)*X.Yd*Yd2Tr*g4[2]
                            + (626./3.)*X.Yd*Yu2Tr*g4[2] + (-204.)*X.Yd*Yd2*z3*g4[2] + (172.)*X.Yd*Yu2*z3*g4[2] + (-24.)*X.Yd*Yd2Tr*z3*g4[2] + (-216.)*X.Yd*Yu2Tr*z3*g4[2] + (833./12.)*X.Yd*g2[0]*g4[2]
                            + (-176./5.)*X.Yd*z3*g2[0]*g4[2] + (987./4.)*X.Yd*g2[1]*g4[2] + (-144.)*X.Yd*z3*g2[1]*g4[2] + (93241./8000.)*X.Yd*g6[0] + (-769./200.)*X.Yd*z3*g6[0] + (455./576.)*X.Yd*g6[1]
                            + (1125./8.)*X.Yd*z3*g6[1] + (-4166./3.)*X.Yd*g6[2] + (640.)*X.Yd*z3*g6[2] + (285./16.)*X.Yd*Yd2*La2 + (-21./16.)*X.Yd*Yu2*La2 + (-135./8.)*X.Yd*Yd2Tr*La2
                            + (-45./8.)*X.Yd*Ye2Tr*La2 + (-135./8.)*X.Yd*Yu2Tr*La2 + (9./4.)*X.Yd*g2[0]*La2 + (45./4.)*X.Yd*g2[1]*La2 + (-9./2.)*X.Yd*La3);
        dX.Ye += loopfactor3*((63./2.)*X.La[0]*X.Ye*Ye4 + (-345./16.)*X.Ye*Ye6 + (45.)*X.La[0]*X.Ye*Ye2*Yd2Tr + (-9./4.)*X.Ye*Ye4*Yd2Tr + (-81./8.)*X.Ye*Ye2*Yd2Tr2 + (45./2.)*X.La[0]*X.Ye*Yd4Tr
                            + (27.)*X.Ye*Ye2*Yd4Tr + (54.)*X.Ye*Yd2Tr*Yd4Tr + (-75./16.)*X.Ye*Yd6Tr + (15.)*X.La[0]*X.Ye*Ye2*Ye2Tr + (-3./4.)*X.Ye*Ye4*Ye2Tr + (-27./4.)*X.Ye*Ye2*Yd2Tr*Ye2Tr
                            + (18.)*X.Ye*Yd4Tr*Ye2Tr + (-9./8.)*X.Ye*Ye2*Ye2Tr2 + (15./2.)*X.La[0]*X.Ye*Ye4Tr + (9.)*X.Ye*Ye2*Ye4Tr + (18.)*X.Ye*Yd2Tr*Ye4Tr + (6.)*X.Ye*Ye2Tr*Ye4Tr + (-25./16.)*X.Ye*Ye6Tr
                            + (39./16.)*X.Ye*Yu2Yd4Tr + (45.)*X.La[0]*X.Ye*Ye2*Yu2Tr + (-9./4.)*X.Ye*Ye4*Yu2Tr + (-81./4.)*X.Ye*Ye2*Yd2Tr*Yu2Tr + (54.)*X.Ye*Yd4Tr*Yu2Tr + (-27./4.)*X.Ye*Ye2*Ye2Tr*Yu2Tr
                            + (18.)*X.Ye*Ye4Tr*Yu2Tr + (-81./8.)*X.Ye*Ye2*Yu2Tr2 + (-9./2.)*X.Ye*Yd2Tr*Yu2Yd2Tr + (-3./2.)*X.Ye*Ye2Tr*Yu2Yd2Tr + (-9./2.)*X.Ye*Yu2Tr*Yu2Yd2Tr + (45./2.)*X.La[0]*X.Ye*Yu4Tr
                            + (27.)*X.Ye*Ye2*Yu4Tr + (54.)*X.Ye*Yd2Tr*Yu4Tr + (18.)*X.Ye*Ye2Tr*Yu4Tr + (54.)*X.Ye*Yu2Tr*Yu4Tr + (39./16.)*X.Ye*Yu4Yd2Tr + (-75./16.)*X.Ye*Yu6Tr + (9./2.)*X.Ye*Ye6*z3
                            + (9.)*X.Ye*Yd6Tr*z3 + (3.)*X.Ye*Ye6Tr*z3 + (9.)*X.Ye*Yu6Tr*z3 + (-99./20.)*X.La[0]*X.Ye*Ye2*g2[0] + (-369./20.)*X.Ye*Ye4*g2[0] + (-33./20.)*X.Ye*Ye2*Yd2Tr*g2[0]
                            + (-81./20.)*X.Ye*Yd2Tr2*g2[0] + (-909./80.)*X.Ye*Yd4Tr*g2[0] + (-171./20.)*X.Ye*Ye2*Ye2Tr*g2[0] + (-27./10.)*X.Ye*Yd2Tr*Ye2Tr*g2[0] + (-9./20.)*X.Ye*Ye2Tr2*g2[0]
                            + (-99./80.)*X.Ye*Ye4Tr*g2[0] + (-177./20.)*X.Ye*Ye2*Yu2Tr*g2[0] + (-81./10.)*X.Ye*Yd2Tr*Yu2Tr*g2[0] + (-27./10.)*X.Ye*Ye2Tr*Yu2Tr*g2[0] + (-81./20.)*X.Ye*Yu2Tr2*g2[0]
                            + (-93./40.)*X.Ye*Yu2Yd2Tr*g2[0] + (-633./80.)*X.Ye*Yu4Tr*g2[0] + (-9./5.)*X.Ye*Ye2*Yd2Tr*z3*g2[0] + (27./5.)*X.Ye*Yd4Tr*z3*g2[0] + (27./5.)*X.Ye*Ye2*Ye2Tr*z3*g2[0]
                            + (-27./5.)*X.Ye*Ye4Tr*z3*g2[0] + (18./5.)*X.Ye*Ye2*Yu2Tr*z3*g2[0] + (24./5.)*X.Ye*Yu2Yd2Tr*z3*g2[0] + (-9./5.)*X.Ye*Yu4Tr*z3*g2[0] + (-135./4.)*X.La[0]*X.Ye*Ye2*g2[1]
                            + (-27./4.)*X.Ye*Ye4*g2[1] + (-81./4.)*X.Ye*Ye2*Yd2Tr*g2[1] + (-81./4.)*X.Ye*Yd2Tr2*g2[1] + (-837./16.)*X.Ye*Yd4Tr*g2[1] + (-27./4.)*X.Ye*Ye2*Ye2Tr*g2[1] + (-27./2.)*X.Ye*Yd2Tr*Ye2Tr*g2[1]
                            + (-9./4.)*X.Ye*Ye2Tr2*g2[1] + (-279./16.)*X.Ye*Ye4Tr*g2[1] + (-81./4.)*X.Ye*Ye2*Yu2Tr*g2[1] + (-81./2.)*X.Ye*Yd2Tr*Yu2Tr*g2[1] + (-27./2.)*X.Ye*Ye2Tr*Yu2Tr*g2[1]
                            + (-81./4.)*X.Ye*Yu2Tr2*g2[1] + (-63./8.)*X.Ye*Yu2Yd2Tr*g2[1] + (-837./16.)*X.Ye*Yu4Tr*g2[1] + (-27.)*X.Ye*Ye2*Yd2Tr*z3*g2[1] + (27.)*X.Ye*Yd4Tr*z3*g2[1] + (-9.)*X.Ye*Ye2*Ye2Tr*z3*g2[1]
                            + (9.)*X.Ye*Ye4Tr*z3*g2[1] + (-27.)*X.Ye*Ye2*Yu2Tr*z3*g2[1] + (27.)*X.Ye*Yu4Tr*z3*g2[1] + (261./80.)*X.La[0]*X.Ye*g2[0]*g2[1] + (-7173./640.)*X.Ye*Ye2*g2[0]*g2[1]
                            + (5469./320.)*X.Ye*Yd2Tr*g2[0]*g2[1] + (2223./320.)*X.Ye*Ye2Tr*g2[0]*g2[1] + (3417./320.)*X.Ye*Yu2Tr*g2[0]*g2[1] + (243./10.)*X.Ye*Ye2*z3*g2[0]*g2[1] + (-27./5.)*X.Ye*Yd2Tr*z3*g2[0]*g2[1]
                            + (54./5.)*X.Ye*Ye2Tr*z3*g2[0]*g2[1] + (-297./10.)*X.Ye*Yu2Tr*z3*g2[0]*g2[1] + (-96.)*X.Ye*Ye2*Yd2Tr*g2[2] + (15./2.)*X.Ye*Yd4Tr*g2[2] + (-96.)*X.Ye*Ye2*Yu2Tr*g2[2]
                            + (57.)*X.Ye*Yu2Yd2Tr*g2[2] + (15./2.)*X.Ye*Yu4Tr*g2[2] + (72.)*X.Ye*Ye2*Yd2Tr*z3*g2[2] + (-72.)*X.Ye*Yd4Tr*z3*g2[2] + (72.)*X.Ye*Ye2*Yu2Tr*z3*g2[2] + (-48.)*X.Ye*Yu2Yd2Tr*z3*g2[2]
                            + (-72.)*X.Ye*Yu4Tr*z3*g2[2] + (-991./60.)*X.Ye*Yd2Tr*g2[0]*g2[2] + (-2419./60.)*X.Ye*Yu2Tr*g2[0]*g2[2] + (12.)*X.Ye*Yd2Tr*z3*g2[0]*g2[2] + (204./5.)*X.Ye*Yu2Tr*z3*g2[0]*g2[2]
                            + (-489./4.)*X.Ye*Yd2Tr*g2[1]*g2[2] + (-489./4.)*X.Ye*Yu2Tr*g2[1]*g2[2] + (108.)*X.Ye*Yd2Tr*z3*g2[1]*g2[2] + (108.)*X.Ye*Yu2Tr*z3*g2[1]*g2[2] + (-621./160.)*X.La[0]*X.Ye*g4[0]
                            + (-105327./6400.)*X.Ye*Ye2*g4[0] + (-54511./9600.)*X.Ye*Yd2Tr*g4[0] + (-31959./3200.)*X.Ye*Ye2Tr*g4[0] + (-33563./1920.)*X.Ye*Yu2Tr*g4[0] + (-1377./200.)*X.Ye*Ye2*z3*g4[0]
                            + (-87./100.)*X.Ye*Yd2Tr*z3*g4[0] + (351./100.)*X.Ye*Ye2Tr*z3*g4[0] + (-3471./100.)*X.Ye*Yu2Tr*z3*g4[0] + (6939./1600.)*X.Ye*g2[1]*g4[0] + (-729./40.)*X.Ye*z3*g2[1]*g4[0]
                            + (7227./100.)*X.Ye*g2[2]*g4[0] + (-396./5.)*X.Ye*z3*g2[2]*g4[0] + (-171./32.)*X.La[0]*X.Ye*g4[1] + (14193./256.)*X.Ye*Ye2*g4[1] + (9099./128.)*X.Ye*Yd2Tr*g4[1]
                            + (3033./128.)*X.Ye*Ye2Tr*g4[1] + (3339./128.)*X.Ye*Yu2Tr*g4[1] + (-135./8.)*X.Ye*Ye2*z3*g4[1] + (-297./4.)*X.Ye*Yd2Tr*z3*g4[1] + (-99./4.)*X.Ye*Ye2Tr*z3*g4[1]
                            + (-243./4.)*X.Ye*Yu2Tr*z3*g4[1] + (2943./320.)*X.Ye*g2[0]*g4[1] + (-243./40.)*X.Ye*z3*g2[0]*g4[1] + (351./4.)*X.Ye*g2[2]*g4[1] + (-108.)*X.Ye*z3*g2[2]*g4[1] + (622./3.)*X.Ye*Yd2Tr*g4[2]
                            + (622./3.)*X.Ye*Yu2Tr*g4[2] + (-24.)*X.Ye*Yd2Tr*z3*g4[2] + (-24.)*X.Ye*Yu2Tr*z3*g4[2] + (607261./8000.)*X.Ye*g6[0] + (-6921./200.)*X.Ye*z3*g6[0] + (455./576.)*X.Ye*g6[1]
                            + (1125./8.)*X.Ye*z3*g6[1] + (285./16.)*X.Ye*Ye2*La2 + (-135./8.)*X.Ye*Yd2Tr*La2 + (-45./8.)*X.Ye*Ye2Tr*La2 + (-135./8.)*X.Ye*Yu2Tr*La2 + (9./4.)*X.Ye*g2[0]*La2
                            + (45./4.)*X.Ye*g2[1]*La2 + (-9./2.)*X.Ye*La3);
        }
        
        if(!X.weylordering || (X.weylordering && X.nloops > 4)) {
        dX.La[0] += loopfactor3*((1440.)*X.La[0]*Yd2Tr*Yd4Tr + (-864.)*Yd4Tr2 + (-5643./4.)*X.La[0]*Yd6Tr + (-891.)*Yd2Tr*Yd6Tr + (156.)*Yd8Tr + (480.)*X.La[0]*Yd4Tr*Ye2Tr + (-297.)*Yd6Tr*Ye2Tr
                               + (480.)*X.La[0]*Yd2Tr*Ye4Tr + (-576.)*Yd4Tr*Ye4Tr + (160.)*X.La[0]*Ye2Tr*Ye4Tr + (-96.)*Ye4Tr2 + (-1881./4.)*X.La[0]*Ye6Tr + (-297.)*Yd2Tr*Ye6Tr + (-99.)*Ye2Tr*Ye6Tr + (52.)*Ye8Tr
                               + (1440.)*X.La[0]*Yd4Tr*Yu2Tr + (-891.)*Yd6Tr*Yu2Tr + (480.)*X.La[0]*Ye4Tr*Yu2Tr + (-297.)*Ye6Tr*Yu2Tr + (126.)*X.La[0]*Yd2Tr*Yu2Yd2Tr + (288.)*Yd4Tr*Yu2Yd2Tr + (42.)*X.La[0]*Ye2Tr*Yu2Yd2Tr
                               + (96.)*Ye4Tr*Yu2Yd2Tr + (126.)*X.La[0]*Yu2Tr*Yu2Yd2Tr + (576.)*Yu2Yd2Tr2 + (132.)*Yu2Yd2Yu2Yd2Tr + (135./4.)*X.La[0]*Yu2Yd4Tr + (135.)*Yd2Tr*Yu2Yd4Tr
                               + (45.)*Ye2Tr*Yu2Yd4Tr + (135.)*Yu2Tr*Yu2Yd4Tr + (-249.)*Yu2Yd6Tr + (1440.)*X.La[0]*Yd2Tr*Yu4Tr + (-1728.)*Yd4Tr*Yu4Tr + (480.)*X.La[0]*Ye2Tr*Yu4Tr + (-576.)*Ye4Tr*Yu4Tr
                               + (1440.)*X.La[0]*Yu2Tr*Yu4Tr + (288.)*Yu2Yd2Tr*Yu4Tr + (-864.)*Yu4Tr2 + (135./4.)*X.La[0]*Yu4Yd2Tr + (135.)*Yd2Tr*Yu4Yd2Tr + (45.)*Ye2Tr*Yu4Yd2Tr + (135.)*Yu2Tr*Yu4Yd2Tr
                               + (750.)*Yu4Yd4Tr + (-5643./4.)*X.La[0]*Yu6Tr +  (-891.)*Yd2Tr*Yu6Tr + (-297.)*Ye2Tr*Yu6Tr + (-891.)*Yu2Tr*Yu6Tr + (-249.)*Yu6Yd2Tr + (156.)*Yu8Tr + (-396.)*X.La[0]*Yd6Tr*z3
                               + (-288.)*Yd8Tr*z3 + (-132.)*X.La[0]*Ye6Tr*z3 + (-96.)*Ye8Tr*z3 + (288.)*X.La[0]*Yu2Yd4Tr*z3 + (-288.)*Yu2Yd6Tr*z3 + (288.)*X.La[0]*Yu4Yd2Tr*z3 + (576.)*Yu4Yd4Tr*z3
                               + (-396.)*X.La[0]*Yu6Tr*z3 + (-288.)*Yu6Yd2Tr*z3 + (-288.)*Yu8Tr*z3 + (-839889./16000.)*g8[0] + (336339./10000.)*z3*g8[0] + (228259./384.)*g8[1] + (-20061./16.)*z3*g8[1]
                               + (-81./5.)*X.La[0]*Yd2Tr2*g2[0] + (-5413./20.)*X.La[0]*Yd4Tr*g2[0] + (5111./20.)*Yd6Tr*g2[0] + (-54./5.)*X.La[0]*Yd2Tr*Ye2Tr*g2[0] + (-9./5.)*X.La[0]*Ye2Tr2*g2[0]
                               + (1557./20.)*X.La[0]*Ye4Tr*g2[0] + (81./4.)*Ye6Tr*g2[0] + (-162./5.)*X.La[0]*Yd2Tr*Yu2Tr*g2[0] + (-54./5.)*X.La[0]*Ye2Tr*Yu2Tr*g2[0] + (-81./5.)*X.La[0]*Yu2Tr2*g2[0]
                               + (-121./2.)*X.La[0]*Yu2Yd2Tr*g2[0] + (-2299./20.)*Yu2Yd4Tr*g2[0] + (-2161./20.)*X.La[0]*Yu4Tr*g2[0] + (1337./20.)*Yu4Yd2Tr*g2[0] + (3467./20.)*Yu6Tr*g2[0]
                               + (1494./5.)*X.La[0]*Yd4Tr*z3*g2[0] + (-120.)*Yd6Tr*z3*g2[0] + (-702./5.)*X.La[0]*Ye4Tr*z3*g2[0] + (792./5.)*Ye6Tr*z3*g2[0] + (-12./5.)*X.La[0]*Yu2Yd2Tr*z3*g2[0]
                               + (624./5.)*Yu2Yd4Tr*z3*g2[0] + (342./5.)*X.La[0]*Yu4Tr*z3*g2[0] + (-672./5.)*Yu4Yd2Tr*z3*g2[0] + (408./5.)*Yu6Tr*z3*g2[0] + (-81.)*X.La[0]*Yd2Tr2*g2[1]
                               + (-4653./4.)*X.La[0]*Yd4Tr*g2[1] + (3411./4.)*Yd6Tr*g2[1] + (-54.)*X.La[0]*Yd2Tr*Ye2Tr*g2[1] + (-9.)*X.La[0]*Ye2Tr2*g2[1] + (-1551./4.)*X.La[0]*Ye4Tr*g2[1] + (1137./4.)*Ye6Tr*g2[1]
                               + (-162.)*X.La[0]*Yd2Tr*Yu2Tr*g2[1] + (-54.)*X.La[0]*Ye2Tr*Yu2Tr*g2[1] + (-81.)*X.La[0]*Yu2Tr2*g2[1] + (-207./2.)*X.La[0]*Yu2Yd2Tr*g2[1] + (477./4.)*Yu2Yd4Tr*g2[1]
                               + (-4653./4.)*X.La[0]*Yu4Tr*g2[1] + (477./4.)*Yu4Yd2Tr*g2[1] + (3411./4.)*Yu6Tr*g2[1] + (1026.)*X.La[0]*Yd4Tr*z3*g2[1] + (-216.)*Yd6Tr*z3*g2[1] + (342.)*X.La[0]*Ye4Tr*z3*g2[1]
                               + (-72.)*Ye6Tr*z3*g2[1] + (108.)*X.La[0]*Yu2Yd2Tr*z3*g2[1] + (1026.)*X.La[0]*Yu4Tr*z3*g2[1] + (-216.)*Yu6Tr*z3*g2[1] + (-9027./80.)*X.La[0]*Yd2Tr*g2[0]*g2[1] + (-81./5.)*Yd2Tr2*g2[0]*g2[1]
                               + (-2591./40.)*Yd4Tr*g2[0]*g2[1] + (-11313./80.)*X.La[0]*Ye2Tr*g2[0]*g2[1] + (-6.)*Yd2Tr*Ye2Tr*g2[0]*g2[1] + (63./5.)*Ye2Tr2*g2[0]*g2[1] + (-549./40.)*Ye4Tr*g2[0]*g2[1]
                               + (-19527./80.)*X.La[0]*Yu2Tr*g2[0]*g2[1] + (-126./5.)*Yd2Tr*Yu2Tr*g2[0]*g2[1] + (174./5.)*Ye2Tr*Yu2Tr*g2[0]*g2[1] + (99./5.)*Yu2Tr2*g2[0]*g2[1] + (301./4.)*Yu2Yd2Tr*g2[0]*g2[1]
                               + (-1871./40.)*Yu4Tr*g2[0]*g2[1] + (72./5.)*X.La[0]*Yd2Tr*z3*g2[0]*g2[1] + (-933./5.)*Yd4Tr*z3*g2[0]*g2[1] + (756./5.)*X.La[0]*Ye2Tr*z3*g2[0]*g2[1] + (-1143./5.)*Ye4Tr*z3*g2[0]*g2[1]
                               + (1062./5.)*X.La[0]*Yu2Tr*z3*g2[0]*g2[1] + (372./5.)*Yu2Yd2Tr*z3*g2[0]*g2[1] + (-2229./5.)*Yu4Tr*z3*g2[0]*g2[1] + (1790.)*X.La[0]*Yd4Tr*g2[2] + (-304.)*Yd6Tr*g2[2]
                               + (164.)*X.La[0]*Yu2Yd2Tr*g2[2] + (-16.)*Yu2Yd4Tr*g2[2] + (1790.)*X.La[0]*Yu4Tr*g2[2] + (-16.)*Yu4Yd2Tr*g2[2] + (-304.)*Yu6Tr*g2[2] + (-2592.)*X.La[0]*Yd4Tr*z3*g2[2] + (1920.)*Yd6Tr*z3*g2[2]
                               + (-192.)*X.La[0]*Yu2Yd2Tr*z3*g2[2] + (-384.)*Yu2Yd4Tr*z3*g2[2] + (-2592.)*X.La[0]*Yu4Tr*z3*g2[2] + (-384.)*Yu4Yd2Tr*z3*g2[2] + (1920.)*Yu6Tr*z3*g2[2] + (-991./15.)*X.La[0]*Yd2Tr*g2[0]*g2[2]
                               + (-2564./15.)*Yd4Tr*g2[0]*g2[2] + (-2419./15.)*X.La[0]*Yu2Tr*g2[0]*g2[2] + (3724./15.)*Yu4Tr*g2[0]*g2[2] + (48.)*X.La[0]*Yd2Tr*z3*g2[0]*g2[2] + (1088./5.)*Yd4Tr*z3*g2[0]*g2[2]
                               + (816./5.)*X.La[0]*Yu2Tr*z3*g2[0]*g2[2] + (-448./5.)*Yu4Tr*z3*g2[0]*g2[2] + (-489.)*X.La[0]*Yd2Tr*g2[1]*g2[2] + (-124.)*Yd4Tr*g2[1]*g2[2] + (-489.)*X.La[0]*Yu2Tr*g2[1]*g2[2]
                               + (-64.)*Yu2Yd2Tr*g2[1]*g2[2] + (-124.)*Yu4Tr*g2[1]*g2[2] + (432.)*X.La[0]*Yd2Tr*z3*g2[1]*g2[2] + (192.)*Yd4Tr*z3*g2[1]*g2[2] + (432.)*X.La[0]*Yu2Tr*z3*g2[1]*g2[2] + (768.)*Yu2Yd2Tr*z3*g2[1]*g2[2]
                               + (192.)*Yu4Tr*z3*g2[1]*g2[2] + (1398./5.)*Yd2Tr*g2[0]*g2[1]*g2[2] + (1494./5.)*Yu2Tr*g2[0]*g2[1]*g2[2] + (-864./5.)*Yd2Tr*z3*g2[0]*g2[1]*g2[2] + (-864./5.)*Yu2Tr*z3*g2[0]*g2[1]*g2[2]
                               + (-149623./2400.)*X.La[0]*Yd2Tr*g4[0] + (-231./50.)*Yd2Tr2*g4[0] + (-98839./1200.)*Yd4Tr*g4[0] + (-44127./800.)*X.La[0]*Ye2Tr*g4[0] + (123./25.)*Yd2Tr*Ye2Tr*g4[0] + (2241./50.)*Ye2Tr2*g4[0]
                               + (10413./80.)*Ye4Tr*g4[0] + (-203887./2400.)*X.La[0]*Yu2Tr*g4[0] + (-51./25.)*Yd2Tr*Yu2Tr*g4[0] + (2103./25.)*Ye2Tr*Yu2Tr*g4[0] + (1857./50.)*Yu2Tr2*g4[0] + (-5973./200.)*Yu2Yd2Tr*g4[0]
                               + (929./48.)*Yu4Tr*g4[0] + (-141./25.)*X.La[0]*Yd2Tr*z3*g4[0] + (-407./10.)*Yd4Tr*z3*g4[0] + (-1107./25.)*X.La[0]*Ye2Tr*z3*g4[0] + (135./2.)*Ye4Tr*z3*g4[0] + (-1347./25.)*X.La[0]*Yu2Tr*z3*g4[0]
                               + (-72./25.)*Yu2Yd2Tr*z3*g4[0] + (2957./50.)*Yu4Tr*z3*g4[0] + (13941./100.)*X.La[0]*g2[1]*g4[0] + (68427./800.)*Yd2Tr*g2[1]*g4[0] + (54153./800.)*Ye2Tr*g2[1]*g4[0] + (76323./800.)*Yu2Tr*g2[1]*g4[0]
                               + (-1323./100.)*X.La[0]*z3*g2[1]*g4[0] + (324./25.)*Yd2Tr*z3*g2[1]*g4[0] + (-108./5.)*Ye2Tr*z3*g2[1]*g4[0] + (-216./25.)*Yu2Tr*z3*g2[1]*g4[0] + (297./5.)*X.La[0]*g2[2]*g4[0]
                               + (2049./25.)*Yd2Tr*g2[2]*g4[0] + (1761./25.)*Yu2Tr*g2[2]*g4[0] + (-1584./25.)*X.La[0]*z3*g2[2]*g4[0] + (-1296./25.)*Yd2Tr*z3*g2[2]*g4[0] + (-1296./25.)*Yu2Tr*z3*g2[2]*g4[0]
                               + (-1683./25.)*g2[1]*g2[2]*g4[0] + (1584./25.)*z3*g2[1]*g2[2]*g4[0] + (-6957./32.)*X.La[0]*Yd2Tr*g4[1] + (27./2.)*Yd2Tr2*g4[1] + (9693./16.)*Yd4Tr*g4[1] + (-2319./32.)*X.La[0]*Ye2Tr*g4[1]
                               + (9.)*Yd2Tr*Ye2Tr*g4[1] + (3./2.)*Ye2Tr2*g4[1] + (3231./16.)*Ye4Tr*g4[1] + (-6957./32.)*X.La[0]*Yu2Tr*g4[1] + (27.)*Yd2Tr*Yu2Tr*g4[1] + (9.)*Ye2Tr*Yu2Tr*g4[1] + (27./2.)*Yu2Tr2*g4[1]
                               + (-2871./8.)*Yu2Yd2Tr*g4[1] + (9693./16.)*Yu4Tr*g4[1] + (-351.)*X.La[0]*Yd2Tr*z3*g4[1] + (-819./2.)*Yd4Tr*z3*g4[1] + (-117.)*X.La[0]*Ye2Tr*z3*g4[1] + (-273./2.)*Ye4Tr*z3*g4[1]
                               + (-351.)*X.La[0]*Yu2Tr*z3*g4[1] + (468.)*Yu2Yd2Tr*z3*g4[1] + (-819./2.)*Yu4Tr*z3*g4[1] + (18411./80.)*X.La[0]*g2[0]*g4[1] + (18297./160.)*Yd2Tr*g2[0]*g4[1] + (4347./160.)*Ye2Tr*g2[0]*g4[1]
                               + (10461./160.)*Yu2Tr*g2[0]*g4[1] + (-1179./20.)*X.La[0]*z3*g2[0]*g4[1] + (216./5.)*Yd2Tr*z3*g2[0]*g4[1] + (-36./5.)*Ye2Tr*z3*g2[0]*g4[1] + (162./5.)*Yu2Tr*z3*g2[0]*g4[1]
                               + (405.)*X.La[0]*g2[2]*g4[1] + (651.)*Yd2Tr*g2[2]*g4[1] + (651.)*Yu2Tr*g2[2]*g4[1] + (-432.)*X.La[0]*z3*g2[2]*g4[1] + (-432.)*Yd2Tr*z3*g2[2]*g4[1] + (-432.)*Yu2Tr*z3*g2[2]*g4[1]
                               + (-459./5.)*g2[0]*g2[2]*g4[1] + (432./5.)*z3*g2[0]*g2[2]*g4[1] + (-81509./1200.)*g4[0]*g4[1] + (19953./200.)*z3*g4[0]*g4[1] + (2488./3.)*X.La[0]*Yd2Tr*g4[2] + (768.)*Yd2Tr2*g4[2]
                               + (-4432./3.)*Yd4Tr*g4[2] + (2488./3.)*X.La[0]*Yu2Tr*g4[2] + (1536.)*Yd2Tr*Yu2Tr*g4[2] + (768.)*Yu2Tr2*g4[2] + (-4432./3.)*Yu4Tr*g4[2] + (-96.)*X.La[0]*Yd2Tr*z3*g4[2] + (256.)*Yd4Tr*z3*g4[2]
                               + (-96.)*X.La[0]*Yu2Tr*z3*g4[2] + (256.)*Yu4Tr*z3*g4[2] + (88639./1000.)*X.La[0]*g6[0] + (145569./4000.)*Yd2Tr*g6[0] + (296163./4000.)*Ye2Tr*g6[0] + (376509./4000.)*Yu2Tr*g6[0]
                               + (-13437./500.)*X.La[0]*z3*g6[0] + (54./25.)*Yd2Tr*z3*g6[0] + (-162./25.)*Ye2Tr*z3*g6[0] + (-108./25.)*Yu2Tr*z3*g6[0] + (-237787./4000.)*g2[1]*g6[0] + (19593./500.)*z3*g2[1]*g6[0]
                               + (-5049./125.)*g2[2]*g6[0] + (4752./125.)*z3*g2[2]*g6[0] + (58031./144.)*X.La[0]*g6[1] + (-6849./32.)*Yd2Tr*g6[1] + (-2283./32.)*Ye2Tr*g6[1] + (-6849./32.)*Yu2Tr*g6[1]
                               + (4419./4.)*X.La[0]*z3*g6[1] + (594.)*Yd2Tr*z3*g6[1] + (198.)*Ye2Tr*z3*g6[1] + (594.)*Yu2Tr*z3*g6[1] + (-33133./144.)*g2[0]*g6[1] + (-243./4.)*z3*g2[0]*g6[1] + (-459.)*g2[2]*g6[1]
                               + (432.)*z3*g2[2]*g6[1] + (-162.)*Yd2Tr2*La2 + (2367./4.)*Yd4Tr*La2 + (-108.)*Yd2Tr*Ye2Tr*La2 + (-18.)*Ye2Tr2*La2 + (789./4.)*Ye4Tr*La2 + (-324.)*Yd2Tr*Yu2Tr*La2
                               + (-108.)*Ye2Tr*Yu2Tr*La2 + (-162.)*Yu2Tr2*La2 + (765./2.)*Yu2Yd2Tr*La2 + (2367./4.)*Yu4Tr*La2 + (378.)*Yd4Tr*z3*La2 + (126.)*Ye4Tr*z3*La2 + (-432.)*Yu2Yd2Tr*z3*La2
                               + (378.)*Yu4Tr*z3*La2 + (1251./40.)*Yd2Tr*g2[0]*La2 + (-1623./40.)*Ye2Tr*g2[0]*La2 + (-117./8.)*Yu2Tr*g2[0]*La2 + (-288./5.)*Yd2Tr*z3*g2[0]*La2 + (144./5.)*Ye2Tr*z3*g2[0]*La2
                               + (-72./5.)*Yu2Tr*z3*g2[0]*La2 + (639./8.)*Yd2Tr*g2[1]*La2 + (213./8.)*Ye2Tr*g2[1]*La2 + (639./8.)*Yu2Tr*g2[1]*La2 + (-216.)*Yd2Tr*z3*g2[1]*La2 + (-72.)*Ye2Tr*z3*g2[1]*La2
                               + (-216.)*Yu2Tr*z3*g2[1]*La2 + (-999./10.)*g2[0]*g2[1]*La2 + (-243./5.)*z3*g2[0]*g2[1]*La2 + (-612.)*Yd2Tr*g2[2]*La2 + (-612.)*Yu2Tr*g2[2]*La2 + (576.)*Yd2Tr*z3*g2[2]*La2
                               + (576.)*Yu2Tr*z3*g2[2]*La2 + (-1881./25.)*g4[0]*La2 + (-729./50.)*z3*g4[0]*La2 + (-1389./16.)*g4[1]*La2 + (-513./2.)*z3*g4[1]*La2 + (873./8.)*Yd2Tr*La3 + (291./8.)*Ye2Tr*La3
                               + (873./8.)*Yu2Tr*La3 + (-237./20.)*g2[0]*La3 + (9./5.)*z3*g2[0]*La3 + (-237./4.)*g2[1]*La3 + (9.)*z3*g2[1]*La3 + (897./8.)*La4 + (63.)*z3*La4);
        }
                               
        if(X.nloops > 3) {
          gauge<3> g9(g7.cwiseProduct(g2));
          // extract quark masses and mixing
          ckm quarks(X.Yu, X.Yd);
          quarks.calculate();
          // extract lepton masses
          yukawa zero;
          pmns leptons(zero, X.Ye);
          leptons.calculate();
          double yt = quarks.get_upyukawas().transpose()[2]; double yt2 = yt*yt; double yt4 = yt2*yt2; double yt6 = yt4*yt2;
          double yb = quarks.get_downyukawas().transpose()[2]; double yb2 = yb*yb; double yb4 = yb2*yb2; double yb6 = yb4*yb2;
          double ytau = leptons.get_elyukawas().transpose()[2]; double ytau2 = ytau*ytau; double ytau4 = ytau2*ytau2; double ytau6 = ytau4*ytau2;

          dX.g[0] += loopfactor4*((39./10.)*La3*g3[0] + (189./20.)*La2*yb2*g3[0] + (-57./20.)*X.La[0]*yb4*g3[0] + (-957./32.)*yb6*g3[0] + (297./20.)*La2*yt2*g3[0] + (-32269./480.)*yb4*yt2*g3[0]
                                + (-237./20.)*X.La[0]*yt4*g3[0] + (-29281./480.)*yb2*yt4*g3[0] + (-13653./160.)*yt6*g3[0] + (183./20.)*La2*ytau2*g3[0] + (-48.)*yb4*ytau2*g3[0] + (-5933./360.)*yb2*yt2*ytau2*g3[0]
                                + (-2371./40.)*yt4*ytau2*g3[0] + (-219./20.)*X.La[0]*ytau4*g3[0] + (-747./20.)*yb2*ytau4*g3[0] + (-1899./40.)*yt2*ytau4*g3[0] + (-4123./160.)*ytau6*g3[0] + (-3./2.)*yb6*z3*g3[0]
                                + (14./5.)*yb4*yt2*z3*g3[0] + (14./5.)*yb2*yt4*z3*g3[0] + (-51./10.)*yt6*z3*g3[0] + (28./15.)*yb2*yt2*ytau2*z3*g3[0] + (-9./2.)*ytau6*z3*g3[0] + (-981./80.)*La2*g2[1]*g3[0]
                                + (-153./20.)*X.La[0]*yb2*g2[1]*g3[0] + (44931./640.)*yb4*g2[1]*g3[0] + (-81./20.)*X.La[0]*yt2*g2[1]*g3[0] + (6131./64.)*yb2*yt2*g2[1]*g3[0] + (71463./640.)*yt4*g2[1]*g3[0]
                                + (-3./20.)*X.La[0]*ytau2*g2[1]*g3[0] + (8607./160.)*yb2*ytau2*g2[1]*g3[0] + (9713./160.)*yt2*ytau2*g2[1]*g3[0] + (31141./640.)*ytau4*g2[1]*g3[0] + (99./20.)*yb4*z3*g2[1]*g3[0]
                                + (381./10.)*yb2*yt2*z3*g2[1]*g3[0] + (639./20.)*yt4*z3*g2[1]*g3[0] + (141./5.)*yb2*ytau2*z3*g2[1]*g3[0] + (57./5.)*yt2*ytau2*z3*g2[1]*g3[0] + (33./4.)*ytau4*z3*g2[1]*g3[0]
                                + (497./20.)*yb4*g2[2]*g3[0] + (6901./90.)*yb2*yt2*g2[2]*g3[0] + (1429./20.)*yt4*g2[2]*g3[0] + (3613./30.)*yb2*ytau2*g2[2]*g3[0] + (3061./30.)*yt2*ytau2*g2[2]*g3[0]
                                + (-12./5.)*yb4*z3*g2[2]*g3[0] + (-1208./15.)*yb2*yt2*z3*g2[2]*g3[0] + (-60.)*yt4*z3*g2[2]*g3[0] + (-336./5.)*yb2*ytau2*z3*g2[2]*g3[0] + (-336./5.)*yt2*ytau2*z3*g2[2]*g3[0]
                                + (-69./5.)*yb2*g2[1]*g2[2]*g3[0] + (367./5.)*yt2*g2[1]*g2[2]*g3[0] + (-138./5.)*yb2*z3*g2[1]*g2[2]*g3[0] + (-474./5.)*yt2*z3*g2[1]*g2[2]*g3[0] + (889./160.)*X.La[0]*g3[0]*g4[1]
                                + (-80311./1280.)*yb2*g3[0]*g4[1] + (-439841./3840.)*yt2*g3[0]*g4[1] + (-305839./3840.)*ytau2*g3[0]*g4[1] + (62./5.)*yb2*z3*g3[0]*g4[1] + (154./5.)*yt2*z3*g3[0]*g4[1]
                                + (121./5.)*ytau2*z3*g3[0]*g4[1] + (-41971./360.)*g2[2]*g3[0]*g4[1] + (1868./15.)*z3*g2[2]*g3[0]*g4[1] + (-83./18.)*yb2*g3[0]*g4[2] + (-5731./90.)*yt2*g3[0]*g4[2]
                                + (-4.)*yb2*z3*g3[0]*g4[2] + (796./5.)*yt2*z3*g3[0]*g4[2] + (-437./3.)*g2[1]*g3[0]*g4[2] + (736./5.)*z3*g2[1]*g3[0]*g4[2] + (-1269./400.)*La2*g5[0] + (-531./100.)*X.La[0]*yb2*g5[0]
                                + (1919./128.)*yb4*g5[0] + (-963./100.)*X.La[0]*yt2*g5[0] + (306401./14400.)*yb2*yt2*g5[0] + (29059./640.)*yt4*g5[0] + (-657./100.)*X.La[0]*ytau2*g5[0] + (11837./2400.)*yb2*ytau2*g5[0]
                                + (70529./2400.)*yt2*ytau2*g5[0] + (89161./3200.)*ytau4*g5[0] + (9./20.)*yb4*z3*g5[0] + (-5./3.)*yb2*yt2*z3*g5[0] + (-357./100.)*yt4*z3*g5[0] + yb2*ytau2*z3*g5[0] + (-274./5.)*yt2*ytau2*z3*g5[0]
                                + (-93./20.)*ytau4*z3*g5[0] + (1917./400.)*X.La[0]*g2[1]*g5[0] + (-7461./640.)*yb2*g2[1]*g5[0] + (-42841./3200.)*yt2*g2[1]*g5[0] + (26109./3200.)*ytau2*g2[1]*g5[0] + (-3./10.)*yb2*z3*g2[1]*g5[0]
                                + (-561./50.)*yt2*z3*g2[1]*g5[0] + (-243./10.)*ytau2*z3*g2[1]*g5[0] + (-79./15.)*yb2*g2[2]*g5[0] + (-503./75.)*yt2*g2[2]*g5[0] + (-18./5.)*yb2*z3*g2[2]*g5[0] + (102./25.)*yt2*z3*g2[2]*g5[0]
                                + (-69./100.)*g2[1]*g2[2]*g5[0] + (572059./57600.)*g4[1]*g5[0] + (-6751./300.)*z3*g4[1]*g5[0] + (83389./675.)*g4[2]*g5[0] + (-68656./225.)*z3*g4[2]*g5[0] + (-117923./11520.)*g3[0]*g6[1]
                                + (-3109./20.)*z3*g3[0]*g6[1] + (1529./15.)*g3[0]*g6[2] + (-4640./9.)*z3*g3[0]*g6[2] + (3627./4000.)*X.La[0]*g7[0] + (3876629./288000.)*yb2*g7[0] + (8978897./288000.)*yt2*g7[0]
                                + (1082567./32000.)*ytau2*g7[0] + (-249./250.)*yb2*z3*g7[0] + (1299./250.)*yt2*z3*g7[0] + (2967./250.)*ytau2*z3*g7[0] + (-3819731./96000.)*g2[1]*g7[0] + (16529./500.)*z3*g2[1]*g7[0]
                                + (-3629273./27000.)*g2[2]*g7[0] + (180076./1125.)*z3*g2[2]*g7[0] + (-143035709./4320000.)*g9[0] + (-1638851./22500.)*z3*g9[0]).real();
          dX.g[1] += loopfactor4*((13./2.)*La3*g3[1] + (75./4.)*La2*yb2*g3[1] + (-39./4.)*X.La[0]*yb4*g3[1] + (-2143./32.)*yb6*g3[1] + (75./4.)*La2*yt2*g3[1] + (-3265./32.)*yb4*yt2*g3[1] + (-39./4.)*X.La[0]*yt4*g3[1]
                                + (-3265./32.)*yb2*yt4*g3[1] + (-2143./32.)*yt6*g3[1] + (25./4.)*La2*ytau2*g3[1] + (-91./4.)*yb4*ytau2*g3[1] + (-113./6.)*yb2*yt2*ytau2*g3[1] + (-91./4.)*yt4*ytau2*g3[1]
                                + (-13./4.)*X.La[0]*ytau4*g3[1] + (-41./2.)*yb2*ytau4*g3[1] + (-41./2.)*yt2*ytau4*g3[1] + (-269./32.)*ytau6*g3[1] + (-9./2.)*yb6*z3*g3[1] + (6.)*yb4*yt2*z3*g3[1] + (6.)*yb2*yt4*z3*g3[1]
                                + (-9./2.)*yt6*z3*g3[1] + (4.)*yb2*yt2*ytau2*z3*g3[1] + (-3./2.)*ytau6*z3*g3[1] + (-327./80.)*La2*g2[0]*g3[1] + (-51./20.)*X.La[0]*yb2*g2[0]*g3[1] + (15937./640.)*yb4*g2[0]*g3[1]
                                + (-27./20.)*X.La[0]*yt2*g2[0]*g3[1] + (33959./960.)*yb2*yt2*g2[0]*g3[1] + (3161./128.)*yt4*g2[0]*g3[1] + (-1./20.)*X.La[0]*ytau2*g2[0]*g3[1] + (5471./480.)*yb2*ytau2*g2[0]*g3[1]
                                + (5401./480.)*yt2*ytau2*g2[0]*g3[1] + (10309./1920.)*ytau4*g2[0]*g3[1] + (99./20.)*yb4*z3*g2[0]*g3[1] + (8.)*yb2*yt2*z3*g2[0]*g3[1] + (153./20.)*yt4*z3*g2[0]*g3[1]
                                + (7./5.)*yb2*ytau2*z3*g2[0]*g3[1] + (-36./5.)*yt2*ytau2*z3*g2[0]*g3[1] + (49./20.)*ytau4*z3*g2[0]*g3[1] + (239./4.)*yb4*g2[2]*g3[1] + (739./6.)*yb2*yt2*g2[2]*g3[1] + (239./4.)*yt4*g2[2]*g3[1]
                                + (33./2.)*yb2*ytau2*g2[2]*g3[1] + (33./2.)*yt2*ytau2*g2[2]*g3[1] + (-36.)*yb4*z3*g2[2]*g3[1] + (-152.)*yb2*yt2*z3*g2[2]*g3[1] + (-36.)*yt4*z3*g2[2]*g3[1] + (-16.)*yb2*ytau2*z3*g2[2]*g3[1]
                                + (-16.)*yt2*ytau2*z3*g2[2]*g3[1] + yb2*g2[0]*g2[2]*g3[1] + (199./15.)*yt2*g2[0]*g2[2]*g3[1] + (-78./5.)*yb2*z3*g2[0]*g2[2]*g3[1] + (-94./5.)*yt2*z3*g2[0]*g2[2]*g3[1]
                                + (457./800.)*X.La[0]*g3[1]*g4[0] + (17903./1280.)*yb2*g3[1]*g4[0] + (465089./19200.)*yt2*g3[1]*g4[0] + (302123./19200.)*ytau2*g3[1]*g4[0] + (-143./50.)*yb2*z3*g3[1]*g4[0]
                                + (-249./50.)*yt2*z3*g3[1]*g4[0] + (-97./50.)*ytau2*z3*g3[1]*g4[0] + (-52297./1800.)*g2[2]*g3[1]*g4[0] + (508./15.)*z3*g2[2]*g3[1]*g4[0] + (-307./6.)*yb2*g3[1]*g4[2] + (-307./6.)*yt2*g3[1]*g4[2]
                                + (84.)*yb2*z3*g3[1]*g4[2] + (84.)*yt2*z3*g3[1]*g4[2] + (-437./9.)*g2[0]*g3[1]*g4[2] + (736./15.)*z3*g2[0]*g3[1]*g4[2] + (-363./16.)*La2*g5[1] + (-75./4.)*X.La[0]*yb2*g5[1]
                                + (30213./128.)*yb4*g5[1] + (-75./4.)*X.La[0]*yt2*g5[1] + (16735./64.)*yb2*yt2*g5[1] + (30213./128.)*yt4*g5[1] + (-25./4.)*X.La[0]*ytau2*g5[1] + (8927./96.)*yb2*ytau2*g5[1]
                                + (7007./96.)*yt2*ytau2*g5[1] + (54931./1152.)*ytau4*g5[1] + (-63./4.)*yb4*z3*g5[1] + (57./2.)*yb2*yt2*z3*g5[1] + (-63./4.)*yt4*z3*g5[1] + (-1.)*yb2*ytau2*z3*g5[1] + (5.)*yt2*ytau2*z3*g5[1]
                                + (-59./12.)*ytau4*z3*g5[1] + (69./16.)*X.La[0]*g2[0]*g5[1] + (-90029./1920.)*yb2*g2[0]*g5[1] + (-102497./1920.)*yt2*g2[0]*g5[1] + (-15341./640.)*ytau2*g2[0]*g5[1] + (-1.)*yb2*z3*g2[0]*g5[1]
                                + (-7.)*yt2*z3*g2[0]*g5[1] + (-24./5.)*ytau2*z3*g2[0]*g5[1] + (-361./3.)*yb2*g2[2]*g5[1] + (-361./3.)*yt2*g2[2]*g5[1] + (-14.)*yb2*z3*g2[2]*g5[1] + (-14.)*yt2*z3*g2[2]*g5[1]
                                + (161./20.)*g2[0]*g2[2]*g5[1] + (-787709./19200.)*g4[0]*g5[1] + (659./100.)*z3*g4[0]*g5[1] + (2587./3.)*g4[2]*g5[1] + (-640.)*z3*g4[2]*g5[1] + (-6418229./288000.)*g3[1]*g6[0]
                                + (21173./1500.)*z3*g3[1]*g6[0] + (257./3.)*g3[1]*g6[2] + (-1760./3.)*z3*g3[1]*g6[2] + (2905./96.)*X.La[0]*g7[1] + (-500665./2304.)*yb2*g7[1] + (-500665./2304.)*yt2*g7[1]
                                + (-500665./6912.)*ytau2*g7[1] + (239./6.)*yb2*z3*g7[1] + (239./6.)*yt2*z3*g7[1] + (239./18.)*ytau2*z3*g7[1] + (-375767./11520.)*g2[0]*g7[1] + (4631./60.)*z3*g2[0]*g7[1]
                                + (-72881./72.)*g2[2]*g7[1] + (4108./3.)*z3*g2[2]*g7[1] + (124660945./62208.)*g9[1] + (-78803./36.)*z3*g9[1]).real();
          dX.g[2] += loopfactor4*((9.)*La2*yb2*g3[2] + (-15.)*X.La[0]*yb4*g3[2] + (-423./4.)*yb6*g3[2] + (9.)*La2*yt2*g3[2] + (-1171./12.)*yb4*yt2*g3[2] + (-15.)*X.La[0]*yt4*g3[2] + (-1171./12.)*yb2*yt4*g3[2]
                                + (-423./4.)*yt6*g3[2] + (-181./8.)*yb4*ytau2*g3[2] + (-515./36.)*yb2*yt2*ytau2*g3[2] + (-181./8.)*yt4*ytau2*g3[2] + (-135./8.)*yb2*ytau4*g3[2] + (-135./8.)*yt2*ytau4*g3[2]
                                + (-6.)*yb6*z3*g3[2] + (-8.)*yb4*yt2*z3*g3[2] + (-8.)*yb2*yt4*z3*g3[2] + (-6.)*yt6*z3*g3[2] + (-16./3.)*yb2*yt2*ytau2*z3*g3[2] + (2869./160.)*yb4*g2[0]*g3[2]
                                + (19033./720.)*yb2*yt2*g2[0]*g3[2] + (3641./160.)*yt4*g2[0]*g3[2] + (4201./240.)*yb2*ytau2*g2[0]*g3[2] + (3649./240.)*yt2*ytau2*g2[0]*g3[2] + (81./10.)*yb4*z3*g2[0]*g3[2]
                                + (-5./3.)*yb2*yt2*z3*g2[0]*g3[2] + (21./10.)*yt4*z3*g2[0]*g3[2] + (-27./5.)*yb2*ytau2*z3*g2[0]*g3[2] + (-27./5.)*yt2*ytau2*z3*g2[0]*g3[2] + (3201./32.)*yb4*g2[1]*g3[2]
                                + (1895./16.)*yb2*yt2*g2[1]*g3[2] + (3201./32.)*yt4*g2[1]*g3[2] + (295./16.)*yb2*ytau2*g2[1]*g3[2] + (295./16.)*yt2*ytau2*g2[1]*g3[2] + (45./2.)*yb4*z3*g2[1]*g3[2]
                                + (3.)*yb2*yt2*z3*g2[1]*g3[2] + (45./2.)*yt4*z3*g2[1]*g3[2] + (9.)*yb2*ytau2*z3*g2[1]*g3[2] + (9.)*yt2*ytau2*z3*g2[1]*g3[2] + (-155./32.)*yb2*g2[0]*g2[1]*g3[2]
                                + (77./160.)*yt2*g2[0]*g2[1]*g3[2] + (-9./2.)*yb2*z3*g2[0]*g2[1]*g3[2] + (-27./2.)*yt2*z3*g2[0]*g2[1]*g3[2] + (210847./14400.)*yb2*g3[2]*g4[0] + (362287./14400.)*yt2*g3[2]*g4[0]
                                + (759./200.)*ytau2*g3[2]*g4[0] + (-7./100.)*yb2*z3*g3[2]*g4[0] + (-19./100.)*yt2*z3*g3[2]*g4[0] + (-46951./4800.)*g2[1]*g3[2]*g4[0] + (973./100.)*z3*g2[1]*g3[2]*g4[0]
                                + (-12887./192.)*yb2*g3[2]*g4[1] + (-12887./192.)*yt2*g3[2]*g4[1] + (63./8.)*ytau2*g3[2]*g4[1] + (117./4.)*yb2*z3*g3[2]*g4[1] + (117./4.)*yt2*z3*g3[2]*g4[1] + (-37597./2880.)*g2[0]*g3[2]*g4[1]
                                + (691./60.)*z3*g2[0]*g3[2]*g4[1] + (427.)*yb4*g5[2] + (4282./9.)*yb2*yt2*g5[2] + (427.)*yt4*g5[2] + (133./3.)*yb2*ytau2*g5[2] + (133./3.)*yt2*ytau2*g5[2] + (-96.)*yb4*z3*g5[2]
                                + (-400./3.)*yb2*yt2*z3*g5[2] + (-96.)*yt4*z3*g5[2] + (-1487./60.)*yb2*g2[0]*g5[2] + (-1283./60.)*yt2*g2[0]*g5[2] + (-104./5.)*yb2*z3*g2[0]*g5[2] + (-8./5.)*yt2*z3*g2[0]*g5[2]
                                + (-473./4.)*yb2*g2[1]*g5[2] + (-473./4.)*yt2*g2[1]*g5[2] + (-72.)*yb2*z3*g2[1]*g5[2] + (-72.)*yt2*z3*g2[1]*g5[2] + (69./20.)*g2[0]*g2[1]*g5[2] + (-17771./270.)*g4[0]*g5[2]
                                + (451./18.)*z3*g4[0]*g5[2] + (953./9.)*g4[1]*g5[2] + (-475./6.)*z3*g4[1]*g5[2] + (-6085099./216000.)*g3[2]*g6[0] + (17473./900.)*z3*g3[2]*g6[0] + (-176815./1728.)*g3[2]*g6[1]
                                + (-935./4.)*z3*g3[2]*g6[1] + (-6709./9.)*yb2*g7[2] + (-6709./9.)*yt2*g7[2] + (272.)*yb2*z3*g7[2] + (272.)*yt2*z3*g7[2] + (-57739./540.)*g2[0]*g7[2] + (8119./45.)*z3*g2[0]*g7[2]
                                + (-5969./12.)*g2[1]*g7[2] + (869.)*z3*g2[1]*g7[2] + (63559./18.)*g9[2] + (-44948./9.)*z3*g9[2]).real();
        }
      }
    }
  } else {
    dX.setZero();
  }
};
