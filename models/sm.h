#include "base.h"

#ifndef SM_H
#define SM_H

class sm : public base<3,1> {
  
 public:
 sm() : base<3,1>() {};
 sm(const gauge<3> g_in, const std::complex<double> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const int nloops_in, const bool weylordering_in) : base<3,1>(g_in, Yu_in, Yd_in, Ye_in, nloops_in, weylordering_in) { La[0] = La_in; };
 sm(const gauge<3> g_in, const std::complex<double> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in, const int nloops_in) : base<3,1>(g_in, Yu_in, Yd_in, Ye_in, nloops_in) { La[0] = La_in; };
 sm(const gauge<3> g_in, const std::complex<double> La_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in) : base<3,1>(g_in, Yu_in, Yd_in, Ye_in) { La[0] = La_in; };
 sm(const base<3,1> &X) : base<3,1>(X) {};
 
  // contains the RGEs
  void operator()(const sm &X, sm &dX, const double);
};

using namespace Eigen;
// define lp infinity norm for sm
namespace boost { namespace numeric { namespace odeint {
      template<>
      struct vector_space_norm_inf< sm > {
	typedef double result_type;
	double operator()( const sm &X ) const {
	  return std::max( std::max( std::max( X.g.lpNorm<Infinity>(), X.Yu.lpNorm<Infinity>()),
				     std::max( X.Yd.lpNorm<Infinity>(), X.Ye.lpNorm<Infinity>())),
			   X.La.lpNorm<Infinity>());
	}
      };
} } }
#endif
