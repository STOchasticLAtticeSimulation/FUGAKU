#include "bias_noise.hpp"
#include <fftw3.h>

std::vector<double> dwlist(double N);
std::vector<double> dwlist_fftw(double N);
std::vector<std::vector<std::vector<std::complex<double>>>> fft_fftw(const std::vector<std::vector<std::vector<std::complex<double>>>>& bk);

int main(int argc, char* argv[]) 
{
  if (argc!=2) {
    std::cout << "Specify the noise file number correctly." << std::endl;
    return -1;
  }

  std::ofstream ofs(noisefilename + std::string(argv[1]) + //std::string(".dat"));
		    std::string(".bin"), std::ios::out | std::ios::binary);//, std::ios::app);
  if (ofs.fail()) {
    std::cout << "The noise file couldn't be opened. 'mkdir noisedata'" << std::endl;
    return -1;
  }

  
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------------

fftw_init_threads();
#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif

  std::cout << "Box size : " << NLnoise << std::endl;
  
  int totalstep = ceil(log((NLnoise/2-1)/sigma)/dN);
  int divstep = int(totalstep/totalnoiseNo);
  int modstep = int(totalstep%totalnoiseNo);
  
  for (int l=0; l<totalnoiseNo+1; l++) {
    std::vector<std::vector<double>> noisedata;
    if (l<totalnoiseNo) {
      noisedata = std::vector<std::vector<double>>(divstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
    } else {
      noisedata = std::vector<std::vector<double>>(modstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
    } 
    int count = 0;
    
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
    for (int i=0; i<divstep; i++) {
      if (l<totalnoiseNo || i<modstep) noisedata[i] = dwlist_fftw((i+l*divstep)*dN);
// #ifdef _OPENMP
// #pragma omp critical
// #endif
      {
        count++;
	std::cout << "\rNoiseGenerating : " << std::setw(3) << 100*count/divstep << "%" << std::flush;
      }
    }
    std::cout << std::endl;

    for (size_t n=0; n<noisedata.size(); n++) {      
      ofs.write(reinterpret_cast<const char*>(noisedata[n].data()), noisedata[n].size()*sizeof(double));
      
      std::cout << "\rExporting :       " << std::setw(3) << 100*n/noisedata.size() << "%" << std::flush;
    }
    std::cout << "\rExporting :       100%" << std::endl;
    
  }

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}


std::vector<double> dwlist(double N) {
  std::vector<std::vector<std::vector<std::complex<double>>>> dwk(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise, 0)));
  int count = 0;
  double nsigma = sigma*exp(N);
  std::vector<double> dwlist(NLnoise*NLnoise*NLnoise,0);
  
  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      if (realpoint(i,j,k,NLnoise)) {
	dwk[i][j][k] = dist(engine);
	count++;
      } else if (complexpoint(i,j,k,NLnoise)) {
	dwk[i][j][k] = (dist(engine) + II*dist(engine))/sqrt(2);
	count++;
      }
    }
  }

  // reflection
  int ip, jp, kp; // reflected index
  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      if (!(realpoint(i,j,k,NLnoise)||complexpoint(i,j,k,NLnoise))) {
	if (i==0) {
	  ip = 0;
	} else {
	  ip = NLnoise-i;
	}
	
	if (j==0) {
	  jp = 0;
	} else {
	  jp = NLnoise-j;
	}
	
	if (k==0) {
	  kp = 0;
	} else {
	  kp = NLnoise-k;
	}
	
	dwk[i][j][k] = conj(dwk[ip][jp][kp]);
	count++;
      }
    }
  }

  if (count==0) {
    return dwlist;
  }
  dwk /= sqrt(count);

  std::vector<std::vector<std::vector<std::complex<double>>>> dwlattice = fft(dwk);
  LOOP{
    dwlist[i*NLnoise*NLnoise + j*NLnoise + k] = dwlattice[i][j][k].real();
  }

  return dwlist;
}



std::vector<double> dwlist_fftw(double N) {
  std::vector<std::vector<std::vector<std::complex<double>>>> dwk(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise, 0)));
  int count = 0;
  double nsigma = sigma*exp(N);
  std::vector<double> dwlist(NLnoise*NLnoise*NLnoise,0);
  
  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      if (realpoint(i,j,k,NLnoise)) {
	dwk[i][j][k] = dist(engine);
	count++;
      } else if (complexpoint(i,j,k,NLnoise)) {
	dwk[i][j][k] = (dist(engine) + II*dist(engine))/sqrt(2);
	count++;
      }
    }
  }

  // reflection
  int ip, jp, kp; // reflected index
  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      if (!(realpoint(i,j,k,NLnoise)||complexpoint(i,j,k,NLnoise))) {
	if (i==0) {
	  ip = 0;
	} else {
	  ip = NLnoise-i;
	}
	
	if (j==0) {
	  jp = 0;
	} else {
	  jp = NLnoise-j;
	}
	
	if (k==0) {
	  kp = 0;
	} else {
	  kp = NLnoise-k;
	}
	
	dwk[i][j][k] = conj(dwk[ip][jp][kp]);
	count++;
      }
    }
  }

  if (count==0) {
    return dwlist;
  }
  dwk /= sqrt(count);

// #ifdef _OPENMP
// #pragma omp critical
// #endif
  {
    std::vector<std::vector<std::vector<std::complex<double>>>> dwlattice = fft_fftw(dwk);
    LOOP{
      dwlist[i*NLnoise*NLnoise + j*NLnoise + k] = dwlattice[i][j][k].real();
    }
  }

  return dwlist;
}



// FFTを行う関数
std::vector<std::vector<std::vector<std::complex<double>>>> fft_fftw(const std::vector<std::vector<std::vector<std::complex<double>>>>& bk) {

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
