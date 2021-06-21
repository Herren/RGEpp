#ifndef NUBASE
#define NUBASE

#include "eigentypes.h"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using namespace Eigen;

template<int n, int m> class nubase {
protected:
  unsigned int nloops;
  
  // use L-(L-1)-(L-2) loop running instead of L-L-L loop running
  bool weylordering;

public:
  // members
  gauge<n> g;         // vector of gauge couplings
  self<m> La;         // Higgs selfcoupling(s)
  yukawa Yu,Yd,Ye,Yn; // Yukawa matrices
  yukawa Ka,Mn;       // Wilson coefficient of the Weinberg operator, neutrino mass matrix       
  
  // constructors
 nubase() : g(), La(), Yu(), Yd(), Ye(), Yn(), Ka(), Mn(), nloops(2), weylordering(false) {};
 nubase(const gauge<n> g_in, const self<m> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const yukawa Yn_in, const yukawa Ka_in, const yukawa Mn_in, const unsigned int nloops_in, const bool weylordering_in)
       : g(g_in), La(La_in), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), Yn(Yn_in), Ka(Ka_in), Mn(Mn_in), nloops(nloops_in), weylordering(weylordering_in) {};
 nubase(const gauge<n> g_in, const self<m> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const yukawa Yn_in, const yukawa Ka_in, const yukawa Mn_in, const unsigned int nloops_in)
       : g(g_in), La(La_in), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), Yn(Yn_in), Ka(Ka_in), Mn(Mn_in), nloops(nloops_in), weylordering(false) {};
 nubase(const gauge<n> g_in, const self<m> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const yukawa Yn_in, const yukawa Ka_in, const yukawa Mn_in)
       : g(g_in), La(La_in), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), Yn(Yn_in), Ka(Ka_in), Mn(Mn_in), nloops(2), weylordering(false) {};
  // constructors that sets all selfcouplings to zero
 nubase(const gauge<n> g_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const yukawa Yn_in, const yukawa Ka_in, const yukawa Mn_in, const unsigned int nloops_in, const bool weylordering_in)
       : g(g_in), La(), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), Yn(Yn_in), Ka(Ka_in), Mn(Mn_in), nloops(nloops_in), weylordering(weylordering_in)  {};
 nubase(const gauge<n> g_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const yukawa Yn_in, const yukawa Ka_in, const yukawa Mn_in, const unsigned int nloops_in)
       : g(g_in), La(), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), Yn(Yn_in), Ka(Ka_in), Mn(Mn_in), nloops(nloops_in), weylordering(false)  {};
 nubase(const gauge<n> g_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const yukawa Yn_in, const yukawa Ka_in, const yukawa Mn_in)
       : g(g_in), La(), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), Yn(Yn_in), Ka(Ka_in), Mn(Mn_in), nloops(2), weylordering(false)  {};
 
  // in-place operations for vector space algebra
  nubase operator+=(const nubase<n,m> &X);
  nubase operator*=(const double a);
  
  // set all parameters to zero:
  void setZero();

  // set the number of loops (default: 2)
  void setNloops(const unsigned int nloops_in);
  
  // enable Weyl ordering (default: false)
  void setWeylordering(const bool weylordering_in);

  // get the number of loops
  unsigned int getNloops() const;
  
  // check if Weyl ordering is enabled
  bool getWeylordering() const;
  
  // check for Landau poles and Nan
  bool check();

  // neutrino member functions

  // extract seesaw scales
  Eigen::Vector3d logthresholds();

  // extract mass matrices in SM and 2HDM
  yukawa getML(double vev);
  yukawa getML(double vev, double tanb);
  
  // integrate out the nth right-handed neutrino
  void integrate_out(int gen);
  
  // macro for eigen member functions
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
};

// in place operations

template<int n, int m> nubase<n,m> nubase<n,m>::operator+=(const nubase<n,m> & X) {
  g += X.g;  La += X.La;  Yu += X.Yu;  Yd += X.Yd;  Ye += X.Ye;  Yn += X.Yn;  Ka += X.Ka;  Mn += X.Mn;
  return *this;
}
  
template<int n, int m> nubase<n,m> nubase<n,m>::operator*=(const double a) {
  g *= a;  La *= a;  Yu *= a;   Yd *= a;   Ye *= a;   Yn *= a;   Ka *= a;   Mn *= a;
  return *this;
}

// set all parameters to zero
template<int n, int m> void nubase<n,m>::setZero() {
  g.setZero();  La.setZero();  Yu.setZero();  Yd.setZero();  Ye.setZero();  Yn.setZero();  Ka.setZero();  Mn.setZero();
}

// set number of loops
template<int n, int m> void nubase<n,m>::setNloops(const unsigned int nloops_in) {
  nloops = nloops_in;
}

// enable Weyl ordering
template<int n, int m> void nubase<n,m>::setWeylordering(const bool weylordering_in) {
  weylordering = weylordering_in;
}

// return number of loops
template<int n, int m> unsigned int nubase<n,m>::getNloops() const {
  return nloops;
}

// return if Weyl ordering is enabled
template<int n, int m> bool nubase<n,m>::getWeylordering() const {
  return weylordering;
}

// check for Landau poles and Nan
template<int n, int m> bool nubase<n,m>::check() {
  return isnotnan(g) && ((&g)->template lpNorm<Infinity>() < 3.5)
    && isnotnan(La) && ((&La)->template lpNorm<Infinity>() < 3.5)
    && isnotnan(Yu) && (Yu.lpNorm<Infinity>() < 3.5)
    && isnotnan(Yd) && (Yd.lpNorm<Infinity>() < 3.5)
    && isnotnan(Ye) && (Ye.lpNorm<Infinity>() < 3.5)
    && isnotnan(Yn) && (Yn.lpNorm<Infinity>() < 3.5)
    && isnotnan(Ka) && (Ka.lpNorm<Infinity>() < 3.5)
    && isnotnan(Mn) && (Mn.lpNorm<Infinity>() < 3.5)
    && (nloops > 0);
}


// Auxiliary functions for decoupling
// extract seesaw scales
template<int n, int m> Vector3d nubase<n,m>::logthresholds() {
  SelfAdjointEigenSolver<yukawa> U(Mn.adjoint()*Mn);
  return U.eigenvalues().array().sqrt().log().matrix();
}

// get mass matrix for left-handed neutrinos
template<int n, int m> yukawa nubase<n,m>::getML(double vev) {
  return -0.5e9*vev*vev*(Ka + 2.*Yn.transpose()*Mn.inverse()*Yn);
}

// get mass matrix for left-handed neutrinos in 2HDM
template<int n, int m> yukawa nubase<n,m>::getML(double vev, double tanb) {
  double vu = sin(atan(tanb))*vev;
  return -0.5e9*vu*vu*(Ka + 2.*Yn.transpose()*Mn.inverse()*Yn);
}

// integrate out the nth right-handed neutrino
template<int n, int m> void nubase<n,m>::integrate_out(int gen) {
  
  // diagonalise Mn
  SelfAdjointEigenSolver<yukawa> U(Mn.adjoint()*Mn);
  Mn = U.eigenvectors().transpose()*Mn*U.eigenvectors();
  //  Mn = U.eigenvalues().asDiagonal();
  
  // rotate Yn
  Yn = U.eigenvectors().transpose()*Yn;
  
  // matching on Weinberg operator
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      Ka(i,j) += 2.*Yn(gen,i)*Yn(gen,j)/Mn(gen,gen);

  // deleting entries in Yn
  Yn.row(gen).setZero();
}


// operations for vector space algebra
template<int n, int m> nubase<n,m> operator+(const nubase<n,m> &lhs, const nubase<n,m> &rhs) {
  return nubase<n,m>(lhs.g+rhs.g, lhs.La+rhs.La, lhs.Yu+rhs.Yu, lhs.Yd+rhs.Yd, lhs.Ye+rhs.Ye, lhs.Yn+rhs.Yn, lhs.Ka+rhs.Ka, lhs.Mn+rhs.Mn, lhs.getNloops(), lhs.getWeylordering());
}

template<int n, int m> nubase<n,m> operator*(const nubase<n,m> &lhs, const double &a) {
  return nubase<n,m>(lhs.g*a, lhs.La*a, lhs.Yu*a, lhs.Yd*a, lhs.Ye*a, lhs.Yn*a , lhs.Ka*a , lhs.Mn*a, lhs.getNloops(), lhs.getWeylordering());
}

template<int n, int m> nubase<n,m> operator*(const double &a, const nubase<n,m> &rhs){
  return nubase<n,m>(a*rhs.g, a*rhs.La, a*rhs.Yu, a*rhs.Yd, a*rhs.Ye, a*rhs.Yn, a*rhs.Ka, a*rhs.Mn, rhs.getNloops(), rhs.getWeylordering());
}

template<int n, int m> nubase<n,m> operator+(const nubase<n,m> &lhs, const double &a) {
  return nubase<n,m>( lhs.g+a, lhs.La+a, lhs.Yu+a, lhs.Yd+a, lhs.Ye+a, lhs.Yn+a, lhs.Ka+a, lhs.Mn+a, lhs.getNloops(), lhs.getWeylordering());
}

template<int n, int m> nubase<n,m> operator+(const double &a, const nubase<n,m> &rhs) {
  return nubase<n,m>(a+rhs.g, a+rhs.La, a+rhs.Yu, a+rhs.Yd, a+rhs.Ye, a+rhs.Yn, a+rhs.Ka, a+rhs.Mn, rhs.getNloops(), rhs.getWeylordering());
}

template<int n, int m> nubase<n,m> operator/(const nubase<n,m> &lhs, const nubase<n,m> &rhs) {
  return nubase<n,m>(lhs.g.cwiseQuotient(rhs.g),
                     lhs.La.cwiseQuotient(rhs.La),
                     lhs.Yu.cwiseQuotient(rhs.Yu),
                     lhs.Yd.cwiseQuotient(rhs.Yd),
		     lhs.Ye.cwiseQuotient(rhs.Ye),
		     lhs.Yn.cwiseQuotient(rhs.Yn),
		     lhs.Ka.cwiseQuotient(rhs.Ka),
		     lhs.Mn.cwiseQuotient(rhs.Mn),
		     lhs.getNloops(),
		     lhs.getWeylordering());
}

template<int n, int m> nubase<n,m> abs(const nubase<n,m> &x) {
  return nubase<n,m>(abs(x.g), abs(x.La), abs(x.Yu), abs(x.Yd), abs(x.Ye), abs(x.Yn), abs(x.Ka), abs(x.Mn), x.getNloops(), x.getWeylordering());
}


// specialisation of vector_space_algebra for Eigen::Matrix
namespace boost { namespace numeric { namespace odeint {
template<typename B,int S1,int S2,int O, int M1, int M2>
	struct algebra_dispatcher< Eigen::Matrix<B,S1,S2,O,M1,M2> >{
	    typedef vector_space_algebra algebra_type;
	};
}}}
#endif
