#ifndef INCLUDED_noise_shell_hpp_
#define INCLUDED_noise_shell_hpp_

#include <cstdint>

// Shared by both FFT backends (src/noise_bias.hpp = FFTW, Mac/default;
// src/noise_bias_ssl2.hpp = Fujitsu C-SSL II, Fugaku/-DUSE_SSL2): which
// k-modes are in the noise shell, which of those are independent degrees
// of freedom vs. Hermitian mirrors, and the per-mode Gaussian draw. None of
// this depends on which library actually runs the FFT.

// |k| (in lattice-site units) for every (nx,ny,nz), precomputed once so
// innsigma() below is a table lookup instead of a fresh sqrt() every call.
inline std::array<double,NLnoiseAll> shellRadius{};

inline void init_shell_radius() {
#pragma clang loop vectorize(disable) interleave(disable)
  for (int nx = 0; nx < NLnoise; nx++) {
    int nxt = (nx<=NLnoise/2 ? nx : nx-NLnoise);
#pragma clang loop vectorize(disable) interleave(disable)
    for (int ny = 0; ny < NLnoise; ny++) {
      int nyt = (ny<=NLnoise/2 ? ny : ny-NLnoise);
#pragma clang loop vectorize(disable) interleave(disable)
      for (int nz = 0; nz < NLnoise; nz++) {
        int nzt = (nz<=NLnoise/2 ? nz : nz-NLnoise);
        shellRadius[nx*NLnoise*NLnoise + ny*NLnoise + nz] = sqrt(nxt*nxt + nyt*nyt + nzt*nzt);
      }
    }
  }
}

// judge if point is in nsigma sphere shell
inline bool innsigma(int nx, int ny, int nz, int Num, double nsigma, double dn) {
  int idx = nx*Num*Num + ny*Num + nz;
  return std::abs(shellRadius[idx] - nsigma) <= dn/2.;
}

// judge real point
inline bool realpoint(int nx, int ny, int nz, int Num) {
  return (nx==0||nx==Num/2) && (ny==0||ny==Num/2) && (nz==0||nz==Num/2);
}

// judge independent complex point
inline bool complexpoint(int nx, int ny, int nz, int Num) {
  int nxt = (nx<=Num/2 ? nx : nx-Num);
  int nyt = (ny<=Num/2 ? ny : ny-Num);
  int nzt = (nz<=Num/2 ? nz : nz-Num);

  return (1<=nxt && nxt!=Num/2 && nyt!=Num/2 && nzt!=Num/2) ||
    (nxt==Num/2 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) || (nxt!=Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt!=Num/2) ||
    (nxt==0 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) ||
    (nxt==Num/2 && nyt==Num/2 && 1<=nzt && nzt!=Num/2) || (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt!=Num/2) ||
    (nxt==0 && 1<=nyt && nyt!=Num/2 && nzt==0) ||
    (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==0) || (1<=nxt && nxt!=Num/2 && nyt==0 && nzt==Num/2) || (nxt==0 && nyt==Num/2 && 1<=nzt && nzt!=Num/2);
}

// Counter-based per-point Gaussian noise for the independent k-modes. Each
// mode's draw(s) depend only on (seed, step, field, its own flat index) --
// no shared mutable RNG state -- so fill loops over this parallelize
// trivially. Changes the exact noise realization vs. the old shared-
// std::mt19937 stream (confirmed acceptable: same seed keeps giving the
// same result, matching pre-existing runs is not required); the shell's
// correlation structure is unaffected since it comes entirely from which
// k-modes get nonzero amplitude, not from how each one's value was drawn.
//
// SplitMix64 (Vigna, public domain): simple, fast, no ~2.5kbit
// std::mt19937 state-init cost per grid point.
inline uint64_t splitmix64_next(uint64_t &state) {
  uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

inline uint64_t hash_mix(uint64_t x) {
  uint64_t state = x;
  return splitmix64_next(state);
}

// Two independent standard normal draws for one lattice point, via
// Box-Muller (2 independent uniforms -> 2 independent normals in one shot
// -- the real+imaginary pair complexpoint needs; realpoint uses z0 only).
inline void point_normals(uint64_t seed, uint64_t step, uint64_t field, uint64_t idx, double &z0, double &z1) {
  uint64_t s = hash_mix(seed);
  s = hash_mix(s ^ step);
  s = hash_mix(s ^ field);
  s = hash_mix(s ^ idx);

  uint64_t r1 = splitmix64_next(s);
  uint64_t r2 = splitmix64_next(s);
  double u1 = ((r1 >> 11) + 1) * (1.0/9007199254740992.0); // (0,1], avoids log(0)
  double u2 = (r2 >> 11) * (1.0/9007199254740992.0);       // [0,1)

  double radius = sqrt(-2.0 * log(u1));
  double theta = 2.0 * M_PI * u2;
  z0 = radius * cos(theta);
  z1 = radius * sin(theta);
}

#endif
