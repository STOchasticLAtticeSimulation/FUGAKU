#ifndef INCLUDED_noise_bias_hpp_
#define INCLUDED_noise_bias_hpp_

// `in` holds the full Hermitian-symmetric spectrum built by dwlist_gen (its
// realpoint/complexpoint selection is not aligned with the axis FFTW halves,
// see init_fftw_global), `inhalf` is the non-redundant r2c/c2r half sliced
// out of it (or filled directly, for biaslist1D's trivially symmetric case),
// and `out` is the real-space result of the c2r transform.
inline fftw_complex *in, *inhalf;
inline double *out;
inline fftw_plan plan;

// |k| (in lattice-site units) for every (nx,ny,nz), precomputed once so
// innsigma() below is a table lookup instead of a fresh sqrt() every call.
// Called from 5 separate NLnoise^3 sweeps per step (dwlist_gen x2, its
// mirror pass, biaslist1D x2), and the sqrt itself doesn't depend on
// nsigma/dn, so it was pure repeated work.
//
// REPRODUCIBILITY: dwlist_gen's first fill loop below is serial specifically
// because its dist(engine) call order (and how many times it's called per
// point -- 0, 1, or 2, depending on innsigma/realpoint/complexpoint) fixes
// the entire noise realization for that step; that loop's compiler can't
// vectorize it anyway (dist(engine) is a stateful, non-vectorizable call
// gating each innsigma check), so it's guaranteed to evaluate innsigma
// scalar, one (i,j,k) at a time, in the exact original iteration order. The
// table below is therefore computed with vectorization explicitly disabled
// too, so every entry is bit-identical to what a fresh inline
// sqrt(nxt*nxt+nyt*nyt+nzt*nzt) would have produced at that call site --
// caching the value can't change which points draw from the RNG, only how
// fast the shell test runs.
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

inline void init_fftw_global() {
  static bool is_initialized = false;
  if (!is_initialized) {
    init_shell_radius();
    fftw_init_threads();
    #ifdef _OPENMP
      fftw_plan_with_nthreads(omp_get_max_threads());
    #endif

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoiseAll);
    inhalf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoiseHalfAll);
    out = fftw_alloc_real(NLnoiseAll);

    if(splan && FFTWwisdom){
      plan = fftw_plan_dft_c2r_3d(NLnoise, NLnoise, NLnoise, inhalf, out, FFTW_PATIENT);
      fftw_export_wisdom_to_filename((sdatadir+"/wisdom"+std::to_string(NLnoise)+"_c2r.dat").c_str());
      std::cout << "Make the FFTW plan." << std::endl;
      FFTwisdomFirst=true;
    }
    else if (FFTWwisdom){
      fftw_import_wisdom_from_filename((sdatadir+"/wisdom"+std::to_string(NLnoise)+"_c2r.dat").c_str());
      plan = fftw_plan_dft_c2r_3d(NLnoise, NLnoise, NLnoise, inhalf, out, FFTW_WISDOM_ONLY);
    }
    else{
      plan = fftw_plan_dft_c2r_3d(NLnoise, NLnoise, NLnoise, inhalf, out, FFTW_MEASURE);
    }

    is_initialized = true;
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
    (nxt==Num/2 && nyt==Num/2 && 1<=nzt && nzt!=Num/2) || (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt==Num/2) ||
    (nxt==0 && 1<=nyt && nyt!=Num/2 && nzt==0) ||
    (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==0) || (1<=nxt && nxt!=Num/2 && nyt==0 && nzt==Num/2) || (nxt==0 && nyt==Num/2 && 1<=nzt && nzt!=Num/2);
}


void dwlist_gen(double N, std::mt19937& engine, int Nfield) {
  int count = 0;
  double nsigma = sigma*exp(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] = 0.0;
    in[i][1] = 0.0;
  }

  for (int i = 0; i < NLnoise; i++) {
    for (int j = 0; j < NLnoise; j++) {
      for (int k = 0; k < NLnoise; k++) {
        int idx = i * NLnoise * NLnoise + j * NLnoise + k;

        if (innsigma(i, j, k, NLnoise, nsigma, dn)) {
          if (realpoint(i, j, k, NLnoise)) {
            in[idx][0] = dist(engine);
            count++;
          } else if (complexpoint(i, j, k, NLnoise)) {
            in[idx][0] = dist(engine) * inv_sqrt2;
            in[idx][1] = dist(engine) * inv_sqrt2;
            count++;
          }
        }
      }
    }
  }


#ifdef _OPENMP
#pragma omp parallel for reduction(+:count)
#endif
    for (int i = 0; i < NLnoise; i++) {
      for (int j = 0; j < NLnoise; j++) {
        for (int k = 0; k < NLnoise; k++) {
          int idx = i * NLnoise * NLnoise + j * NLnoise + k;
          if (innsigma(i, j, k, NLnoise, nsigma, dn)) {
            if (!(realpoint(i, j, k, NLnoise) || complexpoint(i, j, k, NLnoise))) {
              int ip = (i == 0 ? 0 : NLnoise - i);
              int jp = (j == 0 ? 0 : NLnoise - j);
              int kp = (k == 0 ? 0 : NLnoise - k);
              in[idx][0] = in[ip * NLnoise * NLnoise + jp * NLnoise + kp][0];
              in[idx][1] = -in[ip * NLnoise * NLnoise + jp * NLnoise + kp][1];
              count++;
            }
          }
        }
      }
    }

  if (count==0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < NLnoiseAll; i++) {
      dwlist[Nfield][i] = out[i];
    }
    return;
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] /= sqrt(count);
    in[i][1] /= sqrt(count);
  }

  // Slice out the non-redundant half (k <= NLnoise/2); the c2r transform
  // reconstructs the rest via Hermitian symmetry.
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int i = 0; i < NLnoise; i++) {
    for (int j = 0; j < NLnoise; j++) {
      for (int k = 0; k < NLnoiseHalf; k++) {
        int idxfull = i * NLnoise * NLnoise + j * NLnoise + k;
        int idxhalf = i * NLnoise * NLnoiseHalf + j * NLnoiseHalf + k;
        inhalf[idxhalf][0] = in[idxfull][0];
        inhalf[idxhalf][1] = in[idxfull][1];
      }
    }
  }

  fftw_execute(plan);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NLnoiseAll; i++) {
    dwlist[Nfield][i] = out[i];
  }
}


void biaslist1D(double N) {
  int count = 0;
  double nsigma = sigma*exp(N);

  // innsigma depends only on |k|, so the shell indicator is already
  // Hermitian-symmetric: fill the non-redundant half directly, no mirroring
  // needed. count still has to run over the full cube to match the original
  // normalization (1/count summed over the whole shell).
#ifdef _OPENMP
#pragma omp parallel for collapse(3) reduction(+:count)
#endif
  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) count++;
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NLnoiseHalfAll; i++) {
    inhalf[i][0] = 0.0;
    inhalf[i][1] = 0.0;
  }

  if (count==0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < NLnoiseAll; i++) {
      biaslist[0][i] = out[i];
    }
    return;
  }

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int i = 0; i < NLnoise; i++) {
    for (int j = 0; j < NLnoise; j++) {
      for (int k = 0; k < NLnoiseHalf; k++) {
        if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
          int idx = i*NLnoise*NLnoiseHalf + j*NLnoiseHalf + k;
          inhalf[idx][0] = 1.0/count;
        }
      }
    }
  }

  fftw_execute(plan);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NLnoiseAll; i++) {
    biaslist[0][i] = out[i];
  }
}

#endif
