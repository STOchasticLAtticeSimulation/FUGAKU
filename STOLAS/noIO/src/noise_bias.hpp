#ifndef INCLUDED_noise_bias_hpp_
#define INCLUDED_noise_bias_hpp_

inline fftw_complex *in, *out;
inline fftw_plan plan;

inline void init_fftw_global() {
  static bool is_initialized = false;
  if (!is_initialized) {
    fftw_init_threads();
    #ifdef _OPENMP
      fftw_plan_with_nthreads(omp_get_max_threads());
    #endif

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoiseAll);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoiseAll);

    plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    is_initialized = true;
  }
}

// judge if point is in nsigma sphere shell
inline bool innsigma(int nx, int ny, int nz, int Num, double nsigma, double dn) {
  int nxt = (nx<=Num/2 ? nx : nx-Num);
  int nyt = (ny<=Num/2 ? ny : ny-Num);
  int nzt = (nz<=Num/2 ? nz : nz-Num);

  return std::abs(sqrt(nxt*nxt + nyt*nyt + nzt*nzt) - nsigma) <= dn/2.;
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

// #ifdef _OPENMP
// #pragma omp parallel for collapse(3) reduction(+:count)
// #endif
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
      dwlist[0][i] = out[i][0];
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

  fftw_execute(plan);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NLnoiseAll; i++) {
    dwlist[Nfield][i] = out[i][0];
  }
}


void biaslist1D(double N) {
  int count = 0;
  double nsigma = sigma*exp(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] = 0.0;
    in[i][1] = 0.0;
  }

#ifdef _OPENMP
#pragma omp parallel for collapse(3) reduction(+:count)
#endif
  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      int idx = i*NLnoise*NLnoise + j*NLnoise + k;
      in[idx][0] = 1.0;
      in[idx][1] = 0.0;
      count++;
    }
  }

  if (count==0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < NLnoiseAll; i++) {
      biaslist[0][i] = out[i][0];
    }
    return;
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] /= count;
  }

  fftw_execute(plan);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NLnoiseAll; i++) {
    biaslist[0][i] = out[i][0];
  }
}

#endif
