#ifndef INCLUDED_noise_bias_ssl2_hpp_
#define INCLUDED_noise_bias_ssl2_hpp_

#include <cstdlib>
#include "noise_shell.hpp"
#include "cssl.h"

// Fujitsu C-SSL II backend for the c2r transform in dwlist_gen/biaslist1D --
// this system's FFTW build has no SVE (or even NEON) codelets for double
// precision, so fftw_execute() there is pure scalar C and ~25x slower than
// the same transform on a NEON-enabled Mac build. C-SSL II ships an
// SVE-optimized library (libssl2mtsve.a, auto-selected by linking -SSL2
// with -mcpu=a64fx already in CXXFLAGS) with c_dm_v3drcf, a mixed-radix
// (2/3/5/7) 3D real<->complex transform.
//
// Fill/mirror logic is identical to src/noise_bias.hpp (the FFTW backend);
// only the transform call and buffer layout differ, since c_dm_v3drcf is
// in-place and wants its fastest dimension padded to 2*(NLnoise/2+1)
// (matching FFTW's own in-place r2c/c2r convention) rather than a separate
// half-size complex array.

// Full-cube working buffer, same role as FFTW backend's `in`.
inline std::array<std::array<double,2>, NLnoiseAll> workFull{};

// c_dm_v3drcf's in-place buffer: as real data, row i,j occupies
// ssl2buf[row*SSL2_K .. row*SSL2_K+NLnoise-1]; reinterpreted as complex
// pairs (post/pre transform), it's SSL2_K/2 = NLnoiseHalf pairs per row.
constexpr int SSL2_K = 2*NLnoiseHalf;
inline std::array<double, (size_t)NLnoise*NLnoise*SSL2_K> ssl2buf{};

inline void init_fftw_global() {
  static bool is_initialized = false;
  if (!is_initialized) {
    init_shell_radius();
    is_initialized = true;
  }
}

inline void ssl2_c2r_transform() {
  int isin = 1, isn = -1, icon = 0;
  c_dm_v3drcf(ssl2buf.data(), SSL2_K, NLnoise, NLnoise, NLnoise, isin, isn, &icon);
  if (icon != 0) {
    std::cerr << "c_dm_v3drcf failed, icon=" << icon << std::endl;
    std::exit(1);
  }
}

void dwlist_gen(double N, int seed, size_t step, int Nfield) {
  int count = 0;
  double nsigma = sigma*exp(N);

#ifdef _OPENMP
#pragma omp parallel for collapse(3) reduction(+:count)
#endif
  for (int i = 0; i < NLnoise; i++) {
    for (int j = 0; j < NLnoise; j++) {
      for (int k = 0; k < NLnoise; k++) {
        int idx = i * NLnoise * NLnoise + j * NLnoise + k;

        if (innsigma(i, j, k, NLnoise, nsigma, dn) && realpoint(i, j, k, NLnoise)) {
          double z0, z1;
          point_normals((uint64_t)seed, (uint64_t)step, (uint64_t)Nfield, (uint64_t)idx, z0, z1);
          workFull[idx][0] = z0;
          workFull[idx][1] = 0.0;
          count++;
        } else if (innsigma(i, j, k, NLnoise, nsigma, dn) && complexpoint(i, j, k, NLnoise)) {
          double z0, z1;
          point_normals((uint64_t)seed, (uint64_t)step, (uint64_t)Nfield, (uint64_t)idx, z0, z1);
          workFull[idx][0] = z0 * inv_sqrt2;
          workFull[idx][1] = z1 * inv_sqrt2;
          count++;
        } else {
          workFull[idx][0] = 0.0;
          workFull[idx][1] = 0.0;
        }
      }
    }
  }

#ifdef _OPENMP
#pragma omp parallel for collapse(3) reduction(+:count)
#endif
  for (int i = 0; i < NLnoise; i++) {
    for (int j = 0; j < NLnoise; j++) {
      for (int k = 0; k < NLnoise; k++) {
        if (innsigma(i, j, k, NLnoise, nsigma, dn) && !(realpoint(i, j, k, NLnoise) || complexpoint(i, j, k, NLnoise))) {
          int idx = i * NLnoise * NLnoise + j * NLnoise + k;
          int ip = (i == 0 ? 0 : NLnoise - i);
          int jp = (j == 0 ? 0 : NLnoise - j);
          int kp = (k == 0 ? 0 : NLnoise - k);
          int pidx = ip * NLnoise * NLnoise + jp * NLnoise + kp;
          workFull[idx][0] = workFull[pidx][0];
          workFull[idx][1] = -workFull[pidx][1];
          count++;
        }
      }
    }
  }

  if (count==0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < NLnoiseAll; i++) {
      dwlist[Nfield][i] = 0.0;
    }
    return;
  }

  double sqrt_count = sqrt((double)count);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int i = 0; i < NLnoise; i++) {
    for (int j = 0; j < NLnoise; j++) {
      size_t rowFull = (size_t)i*NLnoise*NLnoise + (size_t)j*NLnoise;
      size_t rowSSL2 = ((size_t)i*NLnoise + j) * SSL2_K;
      for (int k = 0; k < NLnoiseHalf; k++) {
        ssl2buf[rowSSL2 + 2*k]     = workFull[rowFull + k][0] / sqrt_count;
        ssl2buf[rowSSL2 + 2*k + 1] = workFull[rowFull + k][1] / sqrt_count;
      }
    }
  }

  ssl2_c2r_transform();

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int i = 0; i < NLnoise; i++) {
    for (int j = 0; j < NLnoise; j++) {
      size_t rowSSL2 = ((size_t)i*NLnoise + j) * SSL2_K;
      size_t rowFlat = (size_t)i*NLnoise*NLnoise + (size_t)j*NLnoise;
      for (int k = 0; k < NLnoise; k++) {
        dwlist[Nfield][rowFlat + k] = ssl2buf[rowSSL2 + k];
      }
    }
  }
}


void biaslist1D(double N) {
  int count = 0;
  double nsigma = sigma*exp(N);

#ifdef _OPENMP
#pragma omp parallel for collapse(3) reduction(+:count)
#endif
  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) count++;
  }

  if (count==0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < NLnoiseAll; i++) {
      biaslist[0][i] = 0.0;
    }
    return;
  }

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int i = 0; i < NLnoise; i++) {
    for (int j = 0; j < NLnoise; j++) {
      size_t rowSSL2 = ((size_t)i*NLnoise + j) * SSL2_K;
      for (int k = 0; k < SSL2_K; k++) ssl2buf[rowSSL2 + k] = 0.0;
      for (int k = 0; k < NLnoiseHalf; k++) {
        if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
          ssl2buf[rowSSL2 + 2*k] = 1.0/count;
        }
      }
    }
  }

  ssl2_c2r_transform();

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int i = 0; i < NLnoise; i++) {
    for (int j = 0; j < NLnoise; j++) {
      size_t rowSSL2 = ((size_t)i*NLnoise + j) * SSL2_K;
      size_t rowFlat = (size_t)i*NLnoise*NLnoise + (size_t)j*NLnoise;
      for (int k = 0; k < NLnoise; k++) {
        biaslist[0][rowFlat + k] = ssl2buf[rowSSL2 + k];
      }
    }
  }
}

#endif
