#include "base.h"

#ifndef THDM_H
#define THDM_H

class thdmi : public base<3,5> {
  
 public:
 thdmi() : base<3,5>() {};
 thdmi(const gauge<3> g_in, const self<5> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const int nloops_in) : base<3,5>(g_in, La_in, Yu_in, Yd_in, Ye_in, nloops_in) {};
 thdmi(const gauge<3> g_in, const self<5> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in) : base<3,5>(g_in, La_in, Yu_in, Yd_in, Ye_in) {};
 thdmi(const base<3,5> &X) : base<3,5>(X) {};
 
  // contains the RGEs
  void operator()(const thdmi &X, thdmi &dX, const double);
};

class thdmii : public base<3,5> {
  
 public:
 thdmii() : base<3,5>() {};
 thdmii(const gauge<3> g_in, const self<5> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const int nloops_in) : base<3,5>(g_in, La_in, Yu_in, Yd_in, Ye_in, nloops_in) {};
 thdmii(const gauge<3> g_in, const self<5> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in) : base<3,5>(g_in, La_in, Yu_in, Yd_in, Ye_in) {};
 thdmii(const base<3,5> &X) : base<3,5>(X) {};
 
  // contains the RGEs
  void operator()(const thdmii &X, thdmii &dX, const double);
};


class thdmx : public base<3,5> {
  
 public:
 thdmx() : base<3,5>() {};
 thdmx(const gauge<3> g_in, const self<5> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const int nloops_in) : base<3,5>(g_in, La_in, Yu_in, Yd_in, Ye_in, nloops_in) {};
 thdmx(const gauge<3> g_in, const self<5> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in) : base<3,5>(g_in, La_in, Yu_in, Yd_in, Ye_in) {};
 thdmx(const base<3,5> &X) : base<3,5>(X) {};
 
  // contains the RGEs
  void operator()(const thdmx &X, thdmx &dX, const double);
};


class thdmy : public base<3,5> {
  
 public:
 thdmy() : base<3,5>() {};
 thdmy(const gauge<3> g_in, const self<5> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const int nloops_in) : base<3,5>(g_in, La_in, Yu_in, Yd_in, Ye_in, nloops_in) {};
 thdmy(const gauge<3> g_in, const self<5> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in) : base<3,5>(g_in, La_in, Yu_in, Yd_in, Ye_in) {};
 thdmy(const base<3,5> &X) : base<3,5>(X) {};
 
  // contains the RGEs
  void operator()(const thdmy &X, thdmy &dX, const double);
};


using namespace Eigen;
// define lp infinity norm for 2HDM type I
namespace boost { namespace numeric { namespace odeint {
      template<>
      struct vector_space_norm_inf< thdmi > {
	typedef double result_type;
	double operator()( const thdmi &X ) const {
	  return std::max( std::max( std::max( X.g.lpNorm<Infinity>(), X.Yu.lpNorm<Infinity>()),
				     std::max( X.Yd.lpNorm<Infinity>(), X.Ye.lpNorm<Infinity>())),
			   X.La.lpNorm<Infinity>());
	}
      };
    } } }

// define lp infinity norm for 2HDM type II
namespace boost { namespace numeric { namespace odeint {
      template<>
      struct vector_space_norm_inf< thdmii > {
	typedef double result_type;
	double operator()( const thdmii &X ) const {
	  return std::max( std::max( std::max( X.g.lpNorm<Infinity>(), X.Yu.lpNorm<Infinity>()),
				     std::max( X.Yd.lpNorm<Infinity>(), X.Ye.lpNorm<Infinity>())),
			   X.La.lpNorm<Infinity>());
	}
      };
    } } }

// define lp infinity norm for 2HDM type X
namespace boost { namespace numeric { namespace odeint {
      template<>
      struct vector_space_norm_inf< thdmx > {
	typedef double result_type;
	double operator()( const thdmx &X ) const {
	  return std::max( std::max( std::max( X.g.lpNorm<Infinity>(), X.Yu.lpNorm<Infinity>()),
				     std::max( X.Yd.lpNorm<Infinity>(), X.Ye.lpNorm<Infinity>())),
			   X.La.lpNorm<Infinity>());
	}
      };
    } } }

// define lp infinity norm for 2HDM type Y
namespace boost { namespace numeric { namespace odeint {
      template<>
      struct vector_space_norm_inf< thdmy > {
	typedef double result_type;
	double operator()( const thdmy &X ) const {
	  return std::max( std::max( std::max( X.g.lpNorm<Infinity>(), X.Yu.lpNorm<Infinity>()),
				     std::max( X.Yd.lpNorm<Infinity>(), X.Ye.lpNorm<Infinity>())),
			   X.La.lpNorm<Infinity>());
	}
      };
    } } }

#endif
