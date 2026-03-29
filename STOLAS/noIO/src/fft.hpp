#ifndef INCLUDED_fft_hpp_
#define INCLUDED_fft_hpp_

#define _USR_MATH_DEFINES
#include <cmath>
#include <complex>
#include <fftw3.h>

std::array<std::complex<double>,NLnoiseAll> bkspectrum{};

void fft_1D_real(const std::array<double,NLnoiseAll>& bk) {
  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);

  for (int i = 0; i < NLnoise; ++i) {
    for (int j = 0; j < NLnoise; ++j) {
      for (int k = 0; k < NLnoise; ++k) {
        int idx = i*NLnoise*NLnoise + j*NLnoise + k;
        in[idx][0] = bk[idx];
        in[idx][1] = 0.;
        idx++;
      }
    }
  }

  fftw_plan plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);


  for (int i = 0; i < NLnoise; ++i) {
    for (int j = 0; j < NLnoise; ++j) {
      for (int k = 0; k < NLnoise; ++k) {
        int idx = i*NLnoise*NLnoise + j*NLnoise + k;
        bkspectrum[idx] = out[idx][0] + II*out[idx][1];
      }
    }
  }

  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);
}


#endif
