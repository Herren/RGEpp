#ifndef BASE
#define BASE

#include "eigentypes.h"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using namespace Eigen;

template<int n, int m> class base {
protected:
  int nloops;
  
  // use L-(L-1)-(L-2) loop running instead of L-L-L loop running
  bool weylordering;

public:
  // members
  gauge<n> g;        // vector of gauge couplings
  self<m> La;        // Higgs selfcoupling(s)
  yukawa<3, 3> Yu,Yd,Ye;   // Yukawa matrices
  
  // constructors
 base() : g(), La(), Yu(), Yd(), Ye(), nloops(2), weylordering(false) {};
 base(const gauge<n> g_in, const self<m> La_in, const yukawa<3,3> Yu_in, const yukawa<3,3> Yd_in, const yukawa<3,3> Ye_in, const int nloops_in, const bool weylordering_in) : g(g_in), La(La_in), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), nloops(nloops_in), weylordering(weylordering_in) {};
 base(const gauge<n> g_in, const self<m> La_in, const yukawa<3,3> Yu_in, const yukawa<3,3> Yd_in, const yukawa<3,3> Ye_in, const int nloops_in) : g(g_in), La(La_in), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), nloops(nloops_in), weylordering(false) {};
 base(const gauge<n> g_in, const self<m> La_in, const yukawa<3,3> Yu_in, const yukawa<3,3> Yd_in, const yukawa<3,3> Ye_in) : g(g_in), La(La_in), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), nloops(2), weylordering(false) {};
  // constructors that sets all selfcouplings to zero
 base(const gauge<n> g_in, const yukawa<3,3> Yu_in, const yukawa<3,3> Yd_in, const yukawa<3,3> Ye_in, const int nloops_in, const bool weylordering_in) : g(g_in), La(), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), nloops(nloops_in), weylordering(weylordering_in) {};
 base(const gauge<n> g_in, const yukawa<3,3> Yu_in, const yukawa<3,3> Yd_in, const yukawa<3,3> Ye_in, const int nloops_in) : g(g_in), La(), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), nloops(nloops_in), weylordering(false) {};
 base(const gauge<n> g_in, const yukawa<3,3> Yu_in, const yukawa<3,3> Yd_in, const yukawa<3,3> Ye_in) : g(g_in), La(), Yu(Yu_in), Yd(Yd_in), Ye(Ye_in), nloops(2), weylordering(false) {};
 
  // in-place operations for vector space algebra
  base operator+=(const base<n,m> &X);
  base operator*=(const double a);
  
  // set all parameters to zero:
  void setZero();

  // set the number of loops (default: 2)
  void setNloops(const int nloops_in);
  
  // enable Weyl ordering (default: false)
  void setWeylordering(const bool weylordering_in);

  // get the number of loops
  unsigned int getNloops() const;

  // check if Weyl ordering is enabled
  bool getWeylordering() const;
  
  // check for Landau poles and Nan
  bool check() const;
  
  // macro for eigen member functions
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
};

// in place operations

template<int n, int m> base<n,m> base<n,m>::operator+=(const base<n,m> & X) {
  g += X.g;  La += X.La;  Yu += X.Yu;  Yd += X.Yd;  Ye += X.Ye;
  return *this;
}
  
template<int n, int m> base<n,m> base<n,m>::operator*=(const double a) {
  g *= a;  La *= a;  Yu *= a;   Yd *= a;   Ye *= a; 
  return *this;
}

// set all parameters to zero
template<int n, int m> void base<n,m>::setZero() {
  g.setZero();  La.setZero();  Yu.setZero();  Yd.setZero();  Ye.setZero();
}

// set number of loops
template<int n, int m> void base<n,m>::setNloops(const int nloops_in) {
  nloops = nloops_in;
}

// enable Weyl ordering
template<int n, int m> void base<n,m>::setWeylordering(const bool weylordering_in) {
  weylordering = weylordering_in;
}

// return number of loops
template<int n, int m> unsigned int base<n,m>::getNloops() const {
  return nloops;
}

// return if Weyl ordering is enabled
template<int n, int m> bool base<n,m>::getWeylordering() const {
  return weylordering;
}


// check for Landau poles and Nan
template<int n, int m> bool base<n,m>::check() const {
  return isnotnan(this->g) && ((&g)->template lpNorm<Infinity>() < 3.5)
    && isnotnan(this->La) && ((&La)->template lpNorm<Infinity>() < 3.5)
    && isnotnan(this->Yu) && ((&Yu)->template lpNorm<Infinity>() < 3.5)
    && isnotnan(this->Yd) && ((&Yd)->template lpNorm<Infinity>() < 3.5)
    && isnotnan(this->Ye) && ((&Ye)->template lpNorm<Infinity>() < 3.5)
    && (nloops > 0);
}


// operations for vector space algebra

template<int n, int m> base<n,m> operator+(const base<n,m> &lhs, const base<n,m> &rhs) {
  return base<n,m>(lhs.g+rhs.g, lhs.La+rhs.La, lhs.Yu+rhs.Yu, lhs.Yd+rhs.Yd, lhs.Ye+rhs.Ye, lhs.getNloops(), lhs.getWeylordering());
}

template<int n, int m> base<n,m> operator*(const base<n,m> &lhs, const double &a) {
  return base<n,m>(lhs.g*a, lhs.La*a, lhs.Yu*a, lhs.Yd*a, lhs.Ye*a, lhs.getNloops(), lhs.getWeylordering());
}

template<int n, int m> base<n,m> operator*(const double &a, const base<n,m> &rhs) {
  return base<n,m>(a*rhs.g, a*rhs.La, a*rhs.Yu, a*rhs.Yd, a*rhs.Ye, rhs.getNloops(), rhs.getWeylordering());
}

template<int n, int m> base<n,m> operator+(const base<n,m> &lhs, const double &a) {
  return base<n,m>(lhs.g+a, lhs.La+a, lhs.Yu+a, lhs.Yd+a, lhs.Ye+a, lhs.getNloops(), lhs.getWeylordering());
}

template<int n, int m> base<n,m> operator+(const double &a , const base<n,m> &rhs) {
  return base<n,m>(a+rhs.g, a+rhs.La, a+rhs.Yu, a+rhs.Yd, a+rhs.Ye, rhs.getNloops(), rhs.getWeylordering());
}

template<int n, int m> base<n,m> operator/(const base<n,m> &lhs, const base<n,m> &rhs) {
  return base<n,m>(lhs.g.cwiseQuotient(rhs.g),
	           lhs.La.cwiseQuotient(rhs.La),
	           lhs.Yu.cwiseQuotient(rhs.Yu),
	           lhs.Yd.cwiseQuotient(rhs.Yd),
	           lhs.Ye.cwiseQuotient(rhs.Ye),
                   lhs.getNloops(),
                   lhs.getWeylordering());
}

template<int n, int m> base<n,m> abs(const base<n,m> &x){
  return base<n,m>(abs(x.g), abs(x.La), abs(x.Yu), abs(x.Yd), abs(x.Ye), x.getNloops(), x.getWeylordering());
}

// specialisation of vector_space_algebra for Eigen::Matrix
namespace boost { namespace numeric { namespace odeint {
template<typename B,int S1,int S2,int O, int M1, int M2>
	struct algebra_dispatcher< Eigen::Matrix<B,S1,S2,O,M1,M2> >{
	    typedef vector_space_algebra algebra_type;
	};
}}}


#endif
