#ifndef EIGENTYPES_H
#define EIGENTYPES_H

#include <complex>
#include <Eigen/Core>

// 3x3 Yukawa matrices
template<int n, int m>
using yukawa = Eigen::Matrix<std::complex<double>, n,m>;

// n-dimensional vector of real gauge couplings
template<int n>
using gauge = Eigen::Matrix<double, n,1>;

// n-dimensional vector of complex selfcouplings
template<int n>
using self = Eigen::Matrix<std::complex<double>, n,1>;

// Component-wise absolute value of Yukawa matrix
template<int n, int m> yukawa<n,m> abs(const yukawa<n,m> &Y) {
  return (Y.cwiseAbs()).template cast<std::complex<double> >();
}

// Check Yukawa matrix for NaN
template<int n, int m> bool isnotnan(const yukawa<n,m> &Y) {
  return Y == Y;
}

// Component-wise absolute value of gauge coupling vector
template<int n> gauge<n> abs(const gauge<n> &g) {
  return g.cwiseAbs();
}

// Component-wise absolute value of selfcoupling vector
template<int n> self<n> abs(const self<n> &La) {
  return (La.cwiseAbs()).template cast<std::complex<double> >();
}

// Check gauge couplings for NaN
template<int n> bool isnotnan(const gauge<n> &g) {
  return g == g;
}

// Check selfcouplings for NaN
template<int n> bool isnotnan(const self<n> &La) {
  return La == La;
}

#endif
