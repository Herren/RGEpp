#include "base.h"

#ifndef MSSM_H
#define MSSM_H

class mssm : public base<3,0> {
  
 public:
 mssm() : base<3,0>() {};
 mssm(const gauge<3> g_in, const yukawa Yu_in, const yukawa Yd_in, const yukawa Ye_in) : base<3,0>(g_in, Yu_in, Yd_in, Ye_in) {};
 mssm(const base<3,0> &X) : base<3,0>(X) {};
 
  // contains the RGEs
  void operator()(const mssm &X, mssm &dX, const double);
};

using namespace Eigen;
// define lp infinity norm for sm
namespace boost { namespace numeric { namespace odeint {
      template<>
      struct vector_space_norm_inf< mssm > {
	typedef double result_type;
	double operator()( const mssm &X ) const {
	  return std::max( std::max( X.g.lpNorm<Infinity>(), X.Yu.lpNorm<Infinity>()),
			   std::max( X.Yd.lpNorm<Infinity>(), X.Ye.lpNorm<Infinity>()));
	}
      };
    } } }

#endif
