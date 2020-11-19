#include "eigentypes.h"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

#ifndef RUNDOWN
#define RUNDOWN

template < typename Stepper, typename Model, typename Observer >
struct rundown{
  rundown(Model & state, double start, double end, double error = 1e-10 , double firststep = 0.01){
    using namespace boost::numeric::odeint;

    // define stepper with error control works only with C++11 standard!
    auto cstepper = make_controlled( error , error , Stepper() );

  // check whether input scale is lower then output scale:
  if (end>start)
    integrate_adaptive(cstepper, Model(), state, start, end, abs(firststep), Observer());

  else{ // start integration across thresholds

  // get thresholds from state type
  Eigen::VectorXd scales= state.logthresholds();

  double here = start; // helping variable

  // integrate to next threshold
  for (int i=scales.size()-1; i>=0; i--){
    if (scales[i] > end){
      integrate_adaptive(cstepper, Model(), state, here, scales[i], -abs(firststep), Observer());
      state.integrate_out(i);
      here = scales[i];
    }
    else integrate_adaptive(cstepper, Model(), state, here, end, -abs(firststep), Observer());
  }
  
  // integrate to end
  if (scales[0] > end)
    integrate_adaptive(cstepper, Model(), state, scales[0], end, -abs(firststep), Observer());
  } // else
  } // contructor
};


// specialisation without observer
template < typename Stepper, typename Model >
struct rundown <Stepper, Model, void>{
  rundown(Model & state, double start, double end, double error = 1e-10 , double firststep = 0.01){
      // define stepper with error control works only with C++11 standard!
    auto cstepper = make_controlled( error , error , Stepper() );
  using namespace boost::numeric::odeint;

  // check whether input scale is lower then output scale:
  if (end>start)
    integrate_adaptive(cstepper, Model(), state, start, end, abs(firststep));

  else{ // start integration across thresholds

  // get thresholds from state type
    Eigen::VectorXd scales= state.logthresholds();

  double here = start; // helping vairable

  // integrate to next threshold
  for (int i=scales.size()-1; i>=0; i--){
    if (scales[i] > end){
      integrate_adaptive(cstepper, Model(), state, here, scales[i], -abs(firststep));
      state.integrate_out(i);
      here = scales[i];
    }
    else integrate_adaptive(cstepper, Model(), state, here, end, -abs(firststep));
  }
  
  // integrate to end
  if (scales[0] > end)
    integrate_adaptive(cstepper, Model(), state, scales[0], end, -abs(firststep));
  } // else
  } // contructor
};


#endif
