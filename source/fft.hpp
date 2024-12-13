#ifndef INCLUDED_fft_hpp_
#define INCLUDED_fft_hpp_

#define _USR_MATH_DEFINES
#include <cmath>
#include <complex>
#include <vector>

#include <fftw3.h>

std::vector<std::complex<double>> dft(std::vector<std::complex<double>> signal); // 1-dim discrete Fourier trs.
std::vector<std::vector<std::complex<double>>> dft(std::vector<std::vector<std::complex<double>>> signal); // 2-dim discrete Fourier trs.
std::vector<std::vector<std::vector<std::complex<double>>>> dft(std::vector<std::vector<std::vector<std::complex<double>>>> signal); // 3-dim discrete Fouriere trs.
std::vector<std::complex<double>> fft(std::vector<std::complex<double>> signal); // 1-dim fast Fourier trs.
std::vector<std::vector<std::complex<double>>> fft(std::vector<std::vector<std::complex<double>>> signal); // 2-dim fast Fourier trs.
std::vector<std::vector<std::vector<std::complex<double>>>> fft(std::vector<std::vector<std::vector<std::complex<double>>>> signal); // 3-dim fast Fourier trs.


// FFTW
std::vector<std::vector<std::vector<std::complex<double>>>> fft_fftw(const std::vector<std::vector<std::vector<std::complex<double>>>>& bk);


// DFT
std::vector<std::complex<double>> dft(std::vector<std::complex<double>> signal) {
  int N = signal.size();
  std::vector<std::complex<double>> spectrum(N);
  
  for (int k = 0; k < N; k++) {
    std::complex<double> sum(0.0, 0.0);
    for (int n = 0; n < N; ++n) {
      double angle = -2*M_PI*k*n/N;
      sum += signal[n] * std::polar(1.0, angle);
    }
    spectrum[k] = sum;
  }
  
  return spectrum;
}

std::vector<std::vector<std::complex<double>>> dft(std::vector<std::vector<std::complex<double>>> signal) {
  // dft index y
  for (std::vector<std::complex<double>> &e : signal) {
    e = dft(e);
  }

  // transpose
  std::vector<std::vector<std::complex<double>>> tmp(signal[0].size(), std::vector<std::complex<double>>(signal.size(),0));
  for (size_t x = 0; x < tmp.size(); x++) {
    for (size_t y = 0; y < tmp[0].size(); y++) {
      tmp[x][y] = signal[y][x];
    }
  }

  // dft index x
  for (std::vector<std::complex<double>> &e : tmp) {
    e = dft(e);
  }

  // transpose
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      signal[x][y] = tmp[y][x];
    }
  }

  return signal;
}

std::vector<std::vector<std::vector<std::complex<double>>>> dft(std::vector<std::vector<std::vector<std::complex<double>>>> signal) {
  // dft indices y and z
  for (std::vector<std::vector<std::complex<double>>> &e : signal) {
    e = dft(e);
  }

  // (x,y,z) to (y,z,x)
  std::vector<std::vector<std::vector<std::complex<double>>>> tmp(signal[0].size(), std::vector<std::vector<std::complex<double>>>(signal[0][0].size(), std::vector<std::complex<double>>(signal.size(),0)));
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      for (size_t z = 0; z < signal[0][0].size(); z++) {
	tmp[y][z][x] = signal[x][y][z];
      }
    }
  }
  
  // dft index x
  for (std::vector<std::vector<std::complex<double>>> &e1 : tmp) {
    for (std::vector<std::complex<double>> &e2 : e1) {
      e2 = dft(e2);
    }
  }

  // (y,z,x) to (x,y,z)
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      for (size_t z = 0; z < signal[0][0].size(); z++) {
	signal[x][y][z] = tmp[y][z][x];
      }
    }
  }

  return signal;
}


// Cooley-Tukey FFT
std::vector<std::complex<double>> fft(std::vector<std::complex<double>> signal) {
  int N = signal.size();
  if (N<2) return signal; // divide until N=2
  
  std::vector<std::complex<double>> even(N/2);
  std::vector<std::complex<double>> odd(N/2);
  for (int i=0; i<N/2; i++) {
    even[i] = signal[2*i];
    odd[i] = signal[2*i+1];
  }
  
  even = fft(even);
  odd = fft(odd);
  
  std::vector<std::complex<double>> spectrum(N);
  for (int k=0; k<N/2; k++) {
    std::complex<double> oddW = std::polar(1.0, -2*M_PI*k/N) * odd[k];
    spectrum[k] = even[k] + oddW;
    spectrum[k+N/2] = even[k] - oddW;
  }
  
  return spectrum;
}

std::vector<std::vector<std::complex<double>>> fft(std::vector<std::vector<std::complex<double>>> signal) {
  // fft index y
  for (std::vector<std::complex<double>> &e : signal) {
    e = fft(e);
  }

  // transpose
  std::vector<std::vector<std::complex<double>>> tmp(signal[0].size(), std::vector<std::complex<double>>(signal.size(),0));
  for (size_t x = 0; x < tmp.size(); x++) {
    for (size_t y = 0; y < tmp[0].size(); y++) {
      tmp[x][y] = signal[y][x];
    }
  }

  // fft index x
  for (std::vector<std::complex<double>> &e : tmp) {
    e = fft(e);
  }

  // transpose
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      signal[x][y] = tmp[y][x];
    }
  }

  return signal;
}

std::vector<std::vector<std::vector<std::complex<double>>>> fft(std::vector<std::vector<std::vector<std::complex<double>>>> signal) {
  // fft indices y and z
  for (std::vector<std::vector<std::complex<double>>> &e : signal) {
    e = fft(e);
  }

  // (x,y,z) to (y,z,x)
  std::vector<std::vector<std::vector<std::complex<double>>>> tmp(signal[0].size(), std::vector<std::vector<std::complex<double>>>(signal[0][0].size(), std::vector<std::complex<double>>(signal.size(),0)));
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      for (size_t z = 0; z < signal[0][0].size(); z++) {
	tmp[y][z][x] = signal[x][y][z];
      }
    }
  }
  
  // fft index x
  for (std::vector<std::vector<std::complex<double>>> &e1 : tmp) {
    for (std::vector<std::complex<double>> &e2 : e1) {
      e2 = fft(e2);
    }
  }

  // (y,z,x) to (x,y,z)
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      for (size_t z = 0; z < signal[0][0].size(); z++) {
	signal[x][y][z] = tmp[y][z][x];
      }
    }
  }

  return signal;
}


// FFTW
std::vector<std::vector<std::vector<std::complex<double>>>> fft_fftw(const std::vector<std::vector<std::vector<std::complex<double>>>>& bk) {

  // FFTW用の入力・出力配列を確保
  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);

  // bk (3D std::vector) -> in (1D fftw_complex) に変換
  int idx = 0;
  for (int i = 0; i < NLnoise; ++i) {
    for (int j = 0; j < NLnoise; ++j) {
      for (int k = 0; k < NLnoise; ++k) {
        in[idx][0] = bk[i][j][k].real();
        in[idx][1] = bk[i][j][k].imag();
        idx++;
      }
    }
  }

  // FFTWプランを作成
  fftw_plan plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  // FFTを実行
  fftw_execute(plan);

  // out (1D fftw_complex) -> biaslattice (3D std::vector) に変換
  std::vector<std::vector<std::vector<std::complex<double>>>> biaslattice(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise)));

  idx = 0;
  for (int i = 0; i < NLnoise; ++i) {
    for (int j = 0; j < NLnoise; ++j) {
      for (int k = 0; k < NLnoise; ++k) {
        biaslattice[i][j][k] = std::complex<double>(out[idx][0], out[idx][1]);
        idx++;
      }
    }
  }


  // メモリ解放
  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);

  return biaslattice;
}


#endif
