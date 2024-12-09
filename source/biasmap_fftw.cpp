#include "bias_noise.hpp"
#include <fftw3.h>

std::vector<double> biaslist(double N);
std::vector<double> biaslist_fftw(double N);
std::vector<std::vector<std::vector<std::complex<double>>>> fft_fftw(const std::vector<std::vector<std::vector<std::complex<double>>>>& bk);

int main() 
{
  std::ofstream ofs(biasfilename);
//   std::ofstream ofs("biasdata/biasmap" + std::string(".bin"), std::ios::out | std::ios::binary);
  if (ofs.fail()) {
    std::cout << "The bias file couldn't be opened. 'mkdir biasdata'" << std::endl;
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
  
  int totalstep = ceil(log((NLnoise/2-1)/sigma)/dN), count = 0;
  // std::vector<std::vector<double>> biasdata(totalstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
  int divstep = int(totalstep/totalnoiseNo);
  int modstep = int(totalstep%totalnoiseNo);
  for (int l=0; l<totalnoiseNo+1; l++) {
    std::vector<std::vector<double>> biasdata;
    std::vector<double> biastemp;
    if (l<totalnoiseNo) {
      biasdata = std::vector<std::vector<double>>(divstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
    } else {
      biasdata = std::vector<std::vector<double>>(modstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
    }
    int count = 0;
  
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
    for (int i=0; i<divstep; i++) {
      if (l<totalnoiseNo || i<modstep) {
        biasdata[i] = biaslist_fftw((i+l*divstep)*dN);
        // biastemp = biaslist_fftw((i+l*divstep)*dN);
        // biasdata[i] = biastemp;
      }
// #ifdef _OPENMP
// #pragma omp critical
// #endif
      {
        // biasdata[i] = biastemp;
        count++;
	std::cout << "\rBiasGenerating : " << std::setw(3) << 100*count/divstep << "%" << std::flush;
      }
    }
    std::cout << std::endl;
    
    for (size_t n=0; n<biasdata.size(); n++) {

      for (size_t i=0; i<biasdata[0].size(); i++) {
        if (l<totalnoiseNo || n<modstep) ofs << biasdata[n][i] << ' ';
      }
      if (l<totalnoiseNo || n<modstep) ofs << std::endl;

    //   ofs.write(reinterpret_cast<const char*>(biasdata[n].data()), biasdata[n].size()*sizeof(double));

      std::cout << "\rExporting :      " << std::setw(3) << 100*n/biasdata.size() << "%" << std::flush;
    }
    std::cout << "\rExporting :      100%" << std::endl;
  }

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}


std::vector<double> biaslist(double N) {
  std::vector<std::vector<std::vector<std::complex<double>>>> bk(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise, 0)));
  int count = 0;
  double nsigma = sigma*exp(N);
  std::vector<double> biaslist(NLnoise*NLnoise*NLnoise,0);
  
  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      bk[i][j][k] = 1;
      count++;
    }
  }

  if (count==0) {
    return biaslist;
  }
  bk /= count;

  std::vector<std::vector<std::vector<std::complex<double>>>> biaslattice = fft(bk);
  LOOP{
    biaslist[i*NLnoise*NLnoise + j*NLnoise + k] = biaslattice[i][j][k].real();
  }

  return biaslist;
}

std::vector<double> biaslist_fftw(double N) {
  std::vector<std::vector<std::vector<std::complex<double>>>> bk(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise, 0)));
  int count = 0;
  double nsigma = sigma*exp(N);
  std::vector<double> biaslist(NLnoise*NLnoise*NLnoise,0);
  
  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      bk[i][j][k] = 1;
      count++;
    }
  }

  if (count==0) {
    return biaslist;
  }
  bk /= count;

  std::vector<std::vector<std::vector<std::complex<double>>>> biaslattice = fft_fftw(bk);
  LOOP{
    biaslist[i*NLnoise*NLnoise + j*NLnoise + k] = biaslattice[i][j][k].real();
  }

  return biaslist;
}


// FFTを行う関数
std::vector<std::vector<std::vector<std::complex<double>>>> fft_fftw(
    const std::vector<std::vector<std::vector<std::complex<double>>>>& bk) {
    
    // データサイズを取得
    // const int NLnoise = bk.size();
    const int N3 = NLnoise * NLnoise * NLnoise;

    // FFTW用の入力・出力配列を確保
    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3);

    // bk (3D std::vector) -> in (1D fftw_complex) に変換
    int idx = 0;
    for (int i = 0; i < NLnoise; ++i) {
        for (int j = 0; j < NLnoise; ++j) {
            for (int k = 0; k < NLnoise; ++k) {
                in[idx][0] = bk[i][j][k].real(); // 実部
                in[idx][1] = bk[i][j][k].imag(); // 虚部
                idx++;
            }
        }
    }

    // FFTWプランを作成
    fftw_plan plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // フーリエ変換を実行
    fftw_execute(plan);

    // out (1D fftw_complex) -> biaslattice (3D std::vector) に変換
    std::vector<std::vector<std::vector<std::complex<double>>>> biaslattice(
        NLnoise, std::vector<std::vector<std::complex<double>>>(
                     NLnoise, std::vector<std::complex<double>>(NLnoise)));

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
    // fftw_cleanup_threads();

    return biaslattice;
}
