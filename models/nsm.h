#include "base.h"

#ifndef NSM_H
#define NSM_H

class nsm : public base<3,3> {
  
 public:
 nsm() : base<3,3>() {};
 nsm(const gauge<3> g_in, std::complex<double> La1_in, std::complex<double> La2_in, std::complex<double> La3_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in)
      : base<3,3>(g_in, Yu_in, Yd_in, Ye_in, 1) { La << La1_in, La2_in, La3_in; };
 nsm(const base<3,3> &X) : base<3,3>(X) {};
 
  // contains the RGEs
  void operator()(const nsm &X, nsm &dX, const double);
};


using namespace Eigen;
// define lp infinity norm for nsm
namespace boost { namespace numeric { namespace odeint {
      template<>
      struct vector_space_norm_inf< nsm > {
	typedef double result_type;
	double operator()( const nsm &X ) const {
	  return std::max( std::max( std::max( X.g.lpNorm<Infinity>(), X.Yu.lpNorm<Infinity>()),
				     std::max( X.Yd.lpNorm<Infinity>(), X.Ye.lpNorm<Infinity>())),
			   X.La.lpNorm<Infinity>());
	}
      };
} } }
#endif
