#include "nubase.h"

#ifndef NUMSSM_H
#define NUMSSM_H

class numssm : public nubase<3,0> {
  
 public:
 numssm() : nubase<3,0>() {};
 numssm(const gauge<3> g_in, yukawa Yu_in, yukawa Yd_in, yukawa Ye_in, yukawa Yn_in, yukawa Ka_in, yukawa Mn_in) : nubase<3,0>(g_in, Yu_in, Yd_in, Ye_in, Yn_in, Ka_in, Mn_in) {};
 numssm(const nubase<3,0> &X) : nubase<3,1>(X) {};

  // RGEs
  void operator() (const numssm & X , numssm & dX, const double);
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
    } } } // namespace

#endif
