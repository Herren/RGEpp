#include "numssm.h"

#define loopfactor 0.006332573977646111 // 1/(4 Pi)^2
#define loopfactor2 0.00004010149318236069 // 1/(4 Pi)^4

// define type for couplings using the vector class
using namespace Eigen;

void numssm::operator() (const numssm & X, numssm & dX, const double) {

  if(check()) {
    // variables for powers of parameters
    gauge<3> g2(X.g.array().square().matrix());
    gauge<3> g3(X.g.array().cube().matrix());
    yukawa Yu2 = X.Yu.adjoint()*X.Yu; std::complex<double> Yu2Tr = Yu2.trace();
    yukawa Yd2 = X.Yd.adjoint()*X.Yd; std::complex<double> Yd2Tr = Yd2.trace();
    yukawa Ye2 = X.Ye.adjoint()*X.Ye; std::complex<double> Ye2Tr = Ye2.trace();
    yukawa Yn2 = X.Yn.adjoint()*X.Yn; std::complex<double> Yn2Tr = Yn2.trace();

    dX.g[0] = ((33./5.)*(g3[0]*loopfactor));
    dX.g[1] = (g3[1]*loopfactor);
    dX.g[2] = (-3.*(g3[2]*loopfactor));
    
    if(!X.weylordering || (X.weylordering && X.nloops > 1)) {
    dX.Yu = (loopfactor*((X.Yu*(((-13./15.)*g2[0])+((-3.*g2[1])+(((-16./3.)*g2[2])+(Yn2Tr+(3.*Yu2Tr))))))+X.Yu*(Yd2+(3.*Yu2))));
    dX.Yd = (loopfactor*((X.Yd*(((-7./15.)*g2[0])+((-3.*g2[1])+(((-16./3.)*g2[2])+((3.*Yd2Tr)+Ye2Tr)))))+X.Yd*((3.*Yd2)+Yu2)));
    dX.Ye = (loopfactor*((X.Ye*(((-9./5.)*g2[0])+((-3.*g2[1])+((3.*Yd2Tr)+Ye2Tr))))+X.Ye*((3.*Ye2)+Yn2)));
    dX.Yn = (loopfactor*((X.Yn*(((-3./5.)*g2[0])+((-3.*g2[1])+(Yn2Tr+(3.*Yu2Tr)))))+X.Yn*(Ye2+(3.*Yn2))));
    dX.Ka = (loopfactor*((X.Ka*(((-6./5.)*g2[0])+((-6.*g2[1])+((2.*Yn2Tr)+(6.*Yu2Tr)))))+(X.Ka*(Ye2+Yn2)+(X.Ye.transpose()*X.Ye.conjugate()+X.Yn.transpose()*X.Yn.conjugate())*X.Ka)));
    dX.Mn = (loopfactor*((2.*X.Mn*X.Yn.conjugate()*X.Yn.transpose())+(2.*X.Yn*X.Yn.adjoint()*X.Mn)));
    }
    
    if(X.nloops > 1) {
      gauge<3> g4(g2.array().square().matrix());
      yukawa Yu4 = Yu2*Yu2;             std::complex<double> Yu4Tr = Yu4.trace();
      yukawa Yd4 = Yd2*Yd2;             std::complex<double> Yd4Tr = Yd4.trace();
      yukawa Ye4 = Ye2*Ye2;             std::complex<double> Ye4Tr = Ye4.trace();
      yukawa Yn4 = Yn2*Yn2;             std::complex<double> Yn4Tr = Yn4.trace();
    
      dX.g[0] += (g3[0]*(loopfactor2*(((199./25.)*g2[0])+(((27./5.)*g2[1])+(((88./5.)*g2[2])+(((-14./5.)*Yd2Tr)+(((-18./5.)*Ye2Tr)+((-26./5.)*Yu2Tr)))))))).real();
      dX.g[1] += (g3[1]*(loopfactor2*(((9./5.)*g2[0])+((25.*g2[1])+((24.*g2[2])+((-6.*Yd2Tr)+((-2.*Ye2Tr)+(-6.*Yu2Tr)))))))).real();
      dX.g[2] += (g3[2]*(loopfactor2*(((11./5.)*g2[0])+((9.*g2[1])+((14.*g2[2])+((-4.*Yd2Tr)+(-4.*Yu2Tr))))))).real();
      
      if(!X.weylordering || (X.weylordering && X.nloops > 2)) {
      dX.Yu += (loopfactor2*(X.Yu*(((2./5.)*(g2[0]*Yd2))+((-3.*(Yd2*Yd2Tr))+((-2.*Yd4)+((-1.*(Yd2*Ye2Tr))+(((2./5.)*(g2[0]*Yu2))+((6.*(g2[1]*Yu2))+((-3.*(Yn2Tr*Yu2))
               +((-9.*(Yu2*Yu2Tr))+((-4.*Yu4)+(-2.*Yd2*Yu2))))))))))+(X.Yu*(((2743./450.)*g4[0])+((g2[0]*g2[1])+(((15./2.)*g4[1])+(((136./45.)*(g2[0]*g2[2]))
               +((8.*(g2[1]*g2[2]))+(((-16./9.)*g4[2])+((-3.*Yn4Tr)+(((4./5.)*(g2[0]*Yu2Tr))+((16.*(g2[2]*Yu2Tr))+((-9.*Yu4Tr)+((-1.*(Ye2*Yn2).trace())
               +(-3.*(X.Yu.adjoint()*X.Yd*X.Yd.adjoint()*X.Yu).trace())))))))))))))));
      dX.Yd += (loopfactor2*(X.Yd*(((4./5.)*(g2[0]*Yd2))+((6.*(g2[1]*Yd2))+((-9.*(Yd2*Yd2Tr))+((-4.*Yd4)+((-3.*(Yd2*Ye2Tr))+(((4./5.)*(g2[0]*Yu2))
               +((-1.*(Yn2Tr*Yu2))+((-3.*(Yu2*Yu2Tr))+((-2.*Yu4)+(-2.*Yu2*Yd2))))))))))+(X.Yd*(((287./90.)*g4[0])+((g2[0]*g2[1])+(((15./2.)*g4[1])
               +(((8./9.)*(g2[0]*g2[2]))+((8.*(g2[1]*g2[2]))+(((-16./9.)*g4[2])+(((-2./5.)*(g2[0]*Yd2Tr))+((16.*(g2[2]*Yd2Tr))+((-9.*Yd4Tr)+(((6./5.)*(g2[0]*Ye2Tr))
               +((-3.*Ye4Tr)+((-1.*(Yn2*Ye2).trace())+(-3.*(X.Yd.adjoint()*X.Yu*X.Yu.adjoint()*X.Yd).trace()))))))))))))))));
      dX.Ye += (loopfactor2*(X.Ye*((6.*(g2[1]*Ye2))+((-9.*(Yd2Tr*Ye2))+((-3.*(Ye2*Ye2Tr))+((-4.*Ye4)+((-1.*(Yn2*Yn2Tr))+((-2.*Yn4)+((-3.*(Yn2*Yu2Tr))+(-2.*Yn2*Ye2))))))))
               +(X.Ye*(((27./2.)*g4[0])+(((9./5.)*(g2[0]*g2[1]))+(((15./2.)*g4[1])+(((-2./5.)*(g2[0]*Yd2Tr))+((16.*(g2[2]*Yd2Tr))+((-9.*Yd4Tr)+(((6./5.)*(g2[0]*Ye2Tr))
               +((-3.*Ye4Tr)+((-1.*(Yn2*Ye2).trace())+(-3.*(X.Yd.adjoint()*X.Yu*X.Yu.adjoint()*X.Yd).trace())))))))))))));
      dX.Yn += (loopfactor2*(X.Yn*(((6./5.)*(g2[0]*Ye2))+((-3.*(Yd2Tr*Ye2))+((-1.*(Ye2*Ye2Tr))+((-2.*Ye4)+(((6./5.)*(g2[0]*Yn2))+((6.*(g2[1]*Yn2))+((-3.*(Yn2*Yn2Tr))
               +((-4.*Yn4)+((-9.*(Yn2*Yu2Tr))+(-2.*Ye2*Yn2))))))))))+(X.Yn*(((207./50.)*g4[0])+(((9./5.)*(g2[0]*g2[1]))+(((15./2.)*g4[1])+((-3.*Yn4Tr)+(((4./5.)*(g2[0]*Yu2Tr))
               +((16.*(g2[2]*Yu2Tr))+((-9.*Yu4Tr)+((-1.*(Ye2*Yn2).trace())+(-3.*(X.Yu.adjoint()*X.Yd*X.Yd.adjoint()*X.Yu).trace()))))))))))));
      dX.Ka += (loopfactor2*(X.Ka*((-2.*Ye4)+((-2.*Yn4)+(((((6./5.)*g2[0])+((-3.*Yd2Tr)+(-1.*Ye2Tr)))*X.Ye.transpose()*X.Ye.conjugate())+(-1.*((Yn2Tr+(3.*Yu2Tr))*X.Yn.transpose()*X.Yn.conjugate())))))
               +((((((6./5.)*g2[0])+((-3.*Yd2Tr)+(-1.*Ye2Tr)))*X.Ye.transpose()*X.Ye.conjugate())+((-1.*((Yn2Tr+(3.*Yu2Tr))*X.Yn.transpose()*X.Yn.conjugate()))
               +((-2.*X.Ye.transpose()*X.Ye.conjugate()*X.Ye.transpose()*X.Ye.conjugate())+(-2.*X.Yn.transpose()*X.Yn.conjugate()*X.Yn.transpose()*X.Yn.conjugate()))))*X.Ka
               +(X.Ka*(((207./25.)*g4[0])+(((18./5.)*(g2[0]*g2[1]))+((15*g4[1])+((-6.*Yn4Tr)+(((8./5.)*(g2[0]*Yu2Tr))+((32.*(g2[2]*Yu2Tr))+((-18.*Yu4Tr)
               +((-2.*(Ye2*Yn2).trace())+(-6.*(X.Yu.adjoint()*X.Yd*X.Yd.adjoint()*X.Yu).trace())))))))))))));
      dX.Mn += (loopfactor2*(X.Mn*(((6./5.)*(g2[0]*X.Yn.conjugate()*X.Yn.transpose()))+((6.*(g2[1]*X.Yn.conjugate()*X.Yn.transpose()))+((-2.*(Yn2Tr*X.Yn.conjugate()*X.Yn.transpose()))
               +((-6.*(Yu2Tr*X.Yn.conjugate()*X.Yn.transpose()))+((-2.*X.Yn.conjugate()*X.Ye.transpose()*X.Ye.conjugate()*X.Yn.transpose())
               +(-2.*X.Yn.conjugate()*X.Yn.transpose()*X.Yn.conjugate()*X.Yn.transpose()))))))+(((6./5.)*(g2[0]*X.Yn*X.Yn.adjoint()))
               +((6.*(g2[1]*X.Yn*X.Yn.adjoint()))+((-2.*(Yn2Tr*X.Yn*X.Yn.adjoint()))+((-6.*(Yu2Tr*X.Yn*X.Yn.adjoint()))+((-2.*X.Yn*Ye2*X.Yn.adjoint())+(-2.*X.Yn*Yn2*X.Yn.adjoint()))))))*X.Mn));
      }
    }
  } else {
    dX.setZero();
  }
};
