#include "nubase.h"

#ifndef NUSM_H
#define NUSM_H

class nusm : public nubase<3,1> {
  
 public:
 nusm() : nubase<3,1>() {};
 nusm(const gauge<3> g_in, std::complex<double> La_in, yukawa Yu_in, yukawa Yd_in, yukawa Ye_in, yukawa Yn_in, yukawa Ka_in, yukawa Mn_in) : nubase<3,1>(g_in, Yu_in, Yd_in, Ye_in, Yn_in, Ka_in, Mn_in) {La[0] = La_in;};
 nusm(const nubase<3,1> &X) : nubase<3,1>(X) {};

  // RGEs
  void operator() (const nusm & X , nusm & dX, const double);
};

using namespace Eigen;
// lp infinity norm for nusm
namespace boost { namespace numeric { namespace odeint {
      template<>
      struct vector_space_norm_inf< nusm >
      {
	typedef double result_type;
	double operator()( const nusm &X ) const
	{
	  return std::max(std::max(
				   std::max( X.g.lpNorm<Infinity>() , X.Yu.lpNorm<Infinity>() ),
				   std::max( X.Yd.lpNorm<Infinity>() , X.Ye.lpNorm<Infinity>() )
				   ),
			  std::max(
				   std::max( X.Yn.lpNorm<Infinity>() , X.Ka.lpNorm<Infinity>() ),
				   std::max( X.Mn.lpNorm<Infinity>() , X.La.lpNorm<Infinity>())
				   )
			  );
	}
      };
    } } } // namespace

#endif
