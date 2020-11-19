#ifndef EIGENTYPES_H
#define EIGENTYPES_H

#include <complex>
#include <Eigen/Core>

using yukawa = Eigen::Matrix<std::complex<double>, 3,3>;

template<int n>
using gauge = Eigen::Matrix<double, n,1>;

template<int n>
using self = Eigen::Matrix<std::complex<double>, n,1>;

inline yukawa abs(const yukawa &Y) {
  return Y.cwiseAbs().cast<std::complex<double> >();
}

inline bool isnotnan(const yukawa &Y) {
  return Y == Y;
}

template<int n> gauge<n> abs(const gauge<n> &g){
  return g.cwiseAbs();
}

template<int n> self<n> abs(const self<n> &La){
  return (La.cwiseAbs()).template cast<std::complex<double> >();
}

template<int n> bool isnotnan(const gauge<n> &g){
  return g == g;
}

template<int n> bool isnotnan(const self<n> &La){
  return La == La;
}

#endif
