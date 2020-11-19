#include "nusm.h"

#define loopfactor 0.006332573977646111 // 1/(4 Pi)^2
#define loopfactor2 0.00004010149318236069 // 1/(4 Pi)^4

// define type for couplings using the vector class
using namespace Eigen;

void nusm::operator() ( const nusm & X , nusm & dX, const double) {

  if (check()){
    // variables for powers of parameters
    gauge<3> g2(X.g.array().square().matrix());
    gauge<3> g3(X.g.array().cube().matrix());
    gauge<3> g4(g2.array().square().matrix());
    std::complex<double> La2 = X.La[0]*X.La[0];

    yukawa Yu2 = X.Yu.adjoint()*X.Yu; std::complex<double> Yu2Tr = Yu2.trace();
    yukawa Yd2 = X.Yd.adjoint()*X.Yd; std::complex<double> Yd2Tr = Yd2.trace();
    yukawa Ye2 = X.Ye.adjoint()*X.Ye; std::complex<double> Ye2Tr = Ye2.trace();
    yukawa Yn2 = X.Yn.adjoint()*X.Yn; std::complex<double> Yn2Tr = Yn2.trace();
    yukawa Yu4 = Yu2*Yu2;             std::complex<double> Yu4Tr = Yu4.trace();
    yukawa Yd4 = Yd2*Yd2;             std::complex<double> Yd4Tr = Yd4.trace();
    yukawa Ye4 = Ye2*Ye2;             std::complex<double> Ye4Tr = Ye4.trace();
    yukawa Yn4 = Yn2*Yn2;             std::complex<double> Yn4Tr = Yn4.trace();

    dX.g[0] = ((41./10.)*(g3[0]*loopfactor));
    dX.g[1] = ((-19./6.)*(g3[1]*loopfactor));
    dX.g[2] = (-7.*(g3[2]*loopfactor));
    dX.Yu = (loopfactor*((X.Yu*(((-17./20.)*g2[0])+(((-9./4.)*g2[1])+((-8.*g2[2])+((3.*Yd2Tr)+(Ye2Tr+(Yn2Tr+(3.*Yu2Tr))))))))+((3./2.)*X.Yu*((-1.*Yd2)+Yu2))));
    dX.Yd = (loopfactor*((X.Yd*(((-1./4.)*g2[0])+(((-9./4.)*g2[1])+((-8.*g2[2])+((3.*Yd2Tr)+(Ye2Tr+(Yn2Tr+(3.*Yu2Tr))))))))+((3./2.)*X.Yd*(Yd2+(-1.*Yu2)))));
    dX.Ye = (loopfactor*((X.Ye*(((-9./4.)*g2[0])+(((-9./4.)*g2[1])+((3.*Yd2Tr)+(Ye2Tr+(Yn2Tr+(3.*Yu2Tr)))))))+((3./2.)*X.Ye*(Ye2+(-1.*Yn2)))));
    dX.Yn = (loopfactor*(((3./2.)*X.Yn*((-1*Ye2)+Yn2))+(((-9./20.)*g2[0])+(((-9./4.)*g2[1])+((3.*Yd2Tr)+(Ye2Tr+(Yn2Tr+(3.*Yu2Tr))))))*X.Yn));
    dX.Ka = (loopfactor*(X.Ka*(((-3./2.)*Ye2)+((1./2.)*Yn2))+(((-3.*g2[1])+(X.La[0]+((6.*Yd2Tr)+((2.*Ye2Tr)+((2.*Yn2Tr)+(6.*Yu2Tr))))))*X.Ka+(((-3./2.)*Ye2.transpose())+((1./2.)*Yn2.transpose()))*X.Ka)));
    dX.Mn = (loopfactor*(X.Mn*Yn2.conjugate()+X.Yn*X.Yn.adjoint()*X.Mn));
    dX.La[0] = (loopfactor*(((27./50.)*g4[0])+(((9./5.)*(g2[0]*g2[1]))+(((9./2.)*g4[1])+(((-9./5.)*(g2[0]*X.La[0]))+((-9.*(g2[1]*X.La[0]))+((6.*La2)+((4.*(X.La[0]*((3.*Yd2Tr)+(Ye2Tr+(Yn2Tr+(3.*Yu2Tr))))))+(-8.*((3.*Yd4Tr)+(Ye4Tr+(Yn4Tr+(3.*Yu4Tr)))))))))))));

    if(X.nloops > 1) {
      gauge<3> g5(g3.cwiseProduct(g2));
      gauge<3> g6(g3.array().square().matrix());
      std::complex<double> La3 = La[0]*La2;
  
      dX.g[0] += (loopfactor2*((199./50.)*g5[0]+(((27./10.)*(g3[0]*g2[1]))+(((44./5.)*(g3[0]*g2[2]))+(((-1./2.)*(g3[0]*abs(Yd2Tr)))+(((-3./2.)*(g3[0]*abs(Ye2Tr)))+((-17./10.)*(g3[0]*abs(Yu2Tr)))))))));
      dX.g[1] += (loopfactor2*((35./6.)*g5[1]+(((9./10.)*(g2[0]*g3[1]))+((12.*(g3[1]*g2[2]))+(((-3./2.)*(g3[1]*abs(Yd2Tr)))+(((-1./2.)*(g3[1]*abs(Ye2Tr)))+((-3./2.)*(g3[1]*abs(Yu2Tr)))))))));
      dX.g[2] += (loopfactor2*(-26.*g5[2]+(((11./10.)*(g2[0]*g3[2]))+(((9./2.)*(g2[1]*g3[2]))+((-2.*(g3[2]*abs(Yd2Tr)))+(-2.*(g3[2]*abs(Yu2Tr))))))));
      dX.Yu += (loopfactor2*(X.Yu*((-1*((((43./80.)*g2[0])+(((-9./16.)*g2[1])+(16.*g2[2])))*Yd2))+(((11./4.)*Yd4)+(((((223./80.)*g2[0])+(((135./16.)*g2[1])+(16.*g2[2])))*Yu2)+((-3.*(X.La[0]*Yu2))+(((3./2.)*Yu4)+(((-1./4.)*Yd2*Yu2)+((-1.*Yu2*Yd2)+((((5./4.)*Yd2)+((-9./4.)*Yu2))*(((3.*Yd2)+(Ye2+(3.*Yu2)))).trace()))))))))+(X.Yu*(((1187./600.)*g4[0])+(((-9./20.)*(g2[0]*g2[1]))+(((-23./4.)*g4[1])+(((19./15.)*(g2[0]*g2[2]))+((9.*(g2[1]*g2[2]))+((-108.*g4[2])+(((3./8.)*La2)+(((5./2.)*(((((1./4.)*g2[0])+(((9./4.)*g2[1])+(8.*g2[2])))*Yd2Tr)+(((3./4.)*((g2[0]+g2[1])*Ye2Tr))+((((17./20.)*g2[0])+(((9./4.)*g2[1])+(8.*g2[2])))*Yu2Tr))))+((-9./4.)*(((3.*Yd4)+(Ye4+((3.*Yu4)+(Yd2*Yu2+((-1./3.)*Yu2*Yd2)))))).trace()))))))))))));
      dX.Yd += (loopfactor2*(X.Yd*(((((187./80.)*g2[0])+(((135./16.)*g2[1])+(16.*g2[2])))*Yd2)+((-3.*(X.La[0]*Yd2))+(((3./2.)*Yd4)+((-1.*((((79./80.)*g2[0])+(((-9./16.)*g2[1])+(16.*g2[2])))*Yu2))+(((11./4.)*Yu4)+((-1.*Yd2*Yu2)+(((-1./4.)*Yu2*Yd2)+((((-9./4.)*Yd2)+((5./4.)*Yu2))*(((3.*Yd2)+(Ye2+(3.*Yu2)))).trace()))))))))+(X.Yd*(((-127./600.)*g4[0])+(((-27./20.)*(g2[0]*g2[1]))+(((-23./4.)*g4[1])+(((31./15.)*(g2[0]*g2[2]))+((9.*(g2[1]*g2[2]))+((-108*g4[2])+(((3./8.)*La2)+(((5./2.)*(((((1./4.)*g2[0])+(((9./4.)*g2[1])+(8.*g2[2])))*Yd2Tr)+(((3./4.)*((g2[0]+g2[1])*Ye2Tr))+((((17./20.)*g2[0])+(((9./4.)*g2[1])+(8.*g2[2])))*Yu2Tr))))+((-9./4.)*(((3.*Yd4)+(Ye4+((3.*Yu4)+(Yd2*Yu2+((-1./3.)*Yu2*Yd2)))))).trace()))))))))))));
      dX.Ye += (loopfactor2*(X.Ye*(((((387./80.)*g2[0])+((135./16.)*g2[1]))*Ye2)+((-3*(X.La[0]*Ye2))+(((3./2.)*Ye4)+((-9./4.)*(Ye2*(((3.*Yd2)+(Ye2+(3.*Yu2)))).trace())))))+(X.Ye*(((1371./200.)*g4[0])+(((-23./4.)*g2[1])+(((27./20.)*(g2[0]*g2[1]))+(((3./8.)*La2)+(((5./2.)*(((((1./4.)*g2[0])+(((9./4.)*g2[1])+(8.*g2[2])))*Yd2Tr)+(((3./4.)*((g2[0]+g2[1])*Ye2Tr))+((((17./20.)*g2[0])+(((9./4.)*g2[1])+(8.*g2[2])))*Yu2Tr))))+((-9./4.)*(((3.*Yd4)+(Ye4+((3.*Yu4)+(Yd2*Yu2+((-1./3.)*Yu2*Yd2)))))).trace())))))))));
      dX.La[0] += (2.*(loopfactor2*(((-3411./1000.)*pow(X.g[0],6))+(((305./8.)*pow(X.g[1],6))+(((-1677./200.)*(g4[0]*g2[1]))+(((-289./40.)*(g2[0]*g4[1]))+(((-1./2.)*((((-1887./200.)*g4[0])+(((-117./20.)*(g2[0]*g2[1]))+((73./8.)*g4[1])))*X.La[0]))+(((1./4.)*((((54./5.)*g2[0])+(54*g2[1]))*La2))+(((-39./4.)*La3)+((g2[0]*(((((9./10.)*g2[0])+((27./5.)*g2[1]))*Yd2Tr)+(((((-9./2.)*g2[0])+((33./5.)*g2[1]))*Ye2Tr)+((((-171./50.)*g2[0])+((63./5.)*g2[1]))*Yu2Tr))))+((5.*(X.La[0]*(((((1./4.)*g2[0])+(((9./4.)*g2[1])+(8.*g2[2])))*Yd2Tr)+(((3./4.)*((g2[0]+g2[1])*Ye2Tr))+((((17./20.)*g2[0])+(((9./4.)*g2[1])+(8.*g2[2])))*Yu2Tr)))))+(((-3./2.)*(g4[1]*(((3.*Yd2)+(Ye2+(3.*Yu2)))).trace()))+((-6.*(La2*(((3.*Yd2)+(Ye2+(3.*Yu2)))).trace()))+((-64.*(g2[2]*((Yd4+Yu4)).trace()))+(((-8./5.)*(g2[0]*(((-1.*Yd4)+((3.*Ye4)+(2.*Yu4)))).trace()))+(((-1./2.)*(X.La[0]*(((3.*Yd4)+(Ye4+(3.*Yu4)))).trace()))+((-21.*(X.La[0]*(Yu2*Yd2).trace()))+((20.*(((3.*Yd4*Yd2)+(Ye4*Ye2+(3.*Yu4*Yu2)))).trace())+(-12.*(Yu2*(Yd2+Yu2)*Yd2).trace())))))))))))))))))));
    }
  } else {
    dX.setZero();
  }
};
