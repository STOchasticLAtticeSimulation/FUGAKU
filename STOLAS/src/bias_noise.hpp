#ifndef INCLUDED_bias_noise_hpp_
#define INCLUDED_bias_noise_hpp_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <sys/time.h>

#include "../parameters.hpp"
#include "fft.hpp"
#include "vec_op.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

// useful macro
#define LOOP for(int i = 0; i < NLnoise; i++) for(int j = 0; j < NLnoise; j++) for(int k = 0; k < NLnoise; k++)
const std::complex<double> II(0,1);

// std::random_device seed;
// std::mt19937 engine(seed());
std::mt19937 engine(seed_val);
std::normal_distribution<> dist(0., 1.);

// judge if point is in nsigma sphere shell
bool innsigma(int nx, int ny, int nz, int Num, double nsigma, double dn) {
  int nxt = (nx<=Num/2 ? nx : nx-Num);
  int nyt = (ny<=Num/2 ? ny : ny-Num);
  int nzt = (nz<=Num/2 ? nz : nz-Num);

  return std::abs(sqrt(nxt*nxt + nyt*nyt + nzt*nzt) - nsigma) <= dn/2.;
}

// judge real point
bool realpoint(int nx, int ny, int nz, int Num) {
  return (nx==0||nx==Num/2) && (ny==0||ny==Num/2) && (nz==0||nz==Num/2);
}

// judge independent complex point
bool complexpoint(int nx, int ny, int nz, int Num) {
  int nxt = (nx<=Num/2 ? nx : nx-Num);
  int nyt = (ny<=Num/2 ? ny : ny-Num);
  int nzt = (nz<=Num/2 ? nz : nz-Num);

  return (1<=nxt && nxt!=Num/2 && nyt!=Num/2 && nzt!=Num/2) ||
    (nxt==Num/2 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) || (nxt!=Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt!=Num/2) ||
    (nxt==0 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) ||
    (nxt==Num/2 && nyt==Num/2 && 1<=nzt && nzt!=Num/2) || (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt==Num/2) ||
    (nxt==0 && 1<=nyt && nyt!=Num/2 && nzt==0) ||
    (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==0) || (1<=nxt && nxt!=Num/2 && nyt==0 && nzt==Num/2) || (nxt==0 && nyt==Num/2 && 1<=nzt && nzt!=Num/2);
}

#endif