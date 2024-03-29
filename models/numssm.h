#include "nubase.h"

#ifndef NUMSSM_H
#define NUMSSM_H

class numssm : public nubase<3,0> {
  
 public:
 numssm() : nubase<3,0>() {};
 numssm(const gauge<3> g_in, const yukawa<3,3> Yu_in, const yukawa<3,3> Yd_in, const yukawa<3,3> Ye_in, const yukawa<3,3> Yn_in, const yukawa<3,3> Ka_in, const yukawa<3,3> Mn_in, const int nloops_in, const bool weylordering_in) : nubase<3,0>(g_in, Yu_in, Yd_in, Ye_in, Yn_in, Ka_in, Mn_in, nloops_in, weylordering_in) {};
 numssm(const gauge<3> g_in, const yukawa<3,3> Yu_in, const yukawa<3,3> Yd_in, const yukawa<3,3> Ye_in, const yukawa<3,3> Yn_in, const yukawa<3,3> Ka_in, const yukawa<3,3> Mn_in) : nubase<3,0>(g_in, Yu_in, Yd_in, Ye_in, Yn_in, Ka_in, Mn_in) {};
 numssm(const nubase<3,0> &X) : nubase<3,0>(X) {};

  // RGEs
  void operator() (const numssm & X, numssm & dX, const double);
};

using namespace Eigen;
// lp infinity norm for numssm
namespace boost { namespace numeric { namespace odeint {
      template<>
      struct vector_space_norm_inf< numssm >
      {
	typedef double result_type;
	double operator()( const numssm &X ) const
	{
	  return std::max(std::max(
				   std::max( X.g.lpNorm<Infinity>() , X.Yu.lpNorm<Infinity>() ),
				   std::max( X.Yd.lpNorm<Infinity>() , X.Ye.lpNorm<Infinity>() )
				   ),
			  std::max(
				   std::max( X.Yn.lpNorm<Infinity>() , X.Ka.lpNorm<Infinity>() ),
				   X.Mn.lpNorm<Infinity>())
			  );
	}
      };
} } }
#endif
