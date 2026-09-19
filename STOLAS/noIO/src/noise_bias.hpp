#ifndef INCLUDED_noise_bias_hpp_
#define INCLUDED_noise_bias_hpp_

#include "noise_shell.hpp"

// `in` holds the full Hermitian-symmetric spectrum built by dwlist_gen (its
// realpoint/complexpoint selection is not aligned with the axis FFTW halves,
// see init_fftw_global), `inhalf` is the non-redundant r2c/c2r half sliced
// out of it (or filled directly, for biaslist1D's trivially symmetric case),
// and `out` is the real-space result of the c2r transform.
inline fftw_complex *in, *inhalf;
inline double *out;
inline fftw_plan plan;

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

void dwlist_gen(double N, int seed, size_t step, int Nfield) {
  int count = 0;
  double nsigma = sigma*exp(N);

  // Fused zero + RNG-fill: independent points (realpoint/complexpoint) draw
  // their own Gaussian(s); everything else (including mirror points, for
  // now) is zeroed. Mirror points are filled in a second pass below by
  // copying their partner's value, NOT by recomputing point_normals() again
  // here -- an earlier version did that (folding all 3 old passes into 1),
  // but on Fugaku it turned out log/cos/sin inside point_normals() is the
  // real cost driver, not the barrier count, so doubling those calls (every
  // mirror point re-running them) canceled out the saved barrier. Keeping
  // this as 2 passes instead of 3 still saves one barrier over the original
  // without doubling the transcendental-function work.
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
          in[idx][0] = z0;
          in[idx][1] = 0.0;
          count++;
        } else if (innsigma(i, j, k, NLnoise, nsigma, dn) && complexpoint(i, j, k, NLnoise)) {
          double z0, z1;
          point_normals((uint64_t)seed, (uint64_t)step, (uint64_t)Nfield, (uint64_t)idx, z0, z1);
          in[idx][0] = z0 * inv_sqrt2;
          in[idx][1] = z1 * inv_sqrt2;
          count++;
        } else {
          in[idx][0] = 0.0;
          in[idx][1] = 0.0;
        }
      }
    }
  }

  // Mirror copy: cheap conjugate copy from the already-written independent
  // partner (guaranteed to be realpoint/complexpoint, never another mirror
  // point -- see innsigma/realpoint/complexpoint above), same as before
  // this function used a pure per-point RNG.
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
          in[idx][0] = in[pidx][0];
          in[idx][1] = -in[pidx][1];
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
      dwlist[Nfield][i] = out[i];
    }
    return;
  }

  // Slice out the non-redundant half (k <= NLnoise/2) while normalizing --
  // the c2r transform reconstructs the rest via Hermitian symmetry, so the
  // redundant half of `in` never needs normalizing at all (was a separate
  // full-cube pass before).
  double sqrt_count = sqrt((double)count);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int i = 0; i < NLnoise; i++) {
    for (int j = 0; j < NLnoise; j++) {
      for (int k = 0; k < NLnoiseHalf; k++) {
        int idxfull = i * NLnoise * NLnoise + j * NLnoise + k;
        int idxhalf = i * NLnoise * NLnoiseHalf + j * NLnoiseHalf + k;
        inhalf[idxhalf][0] = in[idxfull][0] / sqrt_count;
        inhalf[idxhalf][1] = in[idxfull][1] / sqrt_count;
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

  if (count==0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < NLnoiseAll; i++) {
      biaslist[0][i] = out[i];
    }
    return;
  }

  // Fused zero + shell-fill (was two separate full-sweep passes): every
  // point in the non-redundant half either gets the shell weight or zero.
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int i = 0; i < NLnoise; i++) {
    for (int j = 0; j < NLnoise; j++) {
      for (int k = 0; k < NLnoiseHalf; k++) {
        int idx = i*NLnoise*NLnoiseHalf + j*NLnoiseHalf + k;
        inhalf[idx][0] = innsigma(i,j,k,NLnoise,nsigma,dn) ? 1.0/count : 0.0;
        inhalf[idx][1] = 0.0;
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
