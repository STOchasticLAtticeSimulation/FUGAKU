#include "bias_noise.hpp"

std::vector<double> biaslist(double N);
std::vector<double> biaslist_fftw(double N);

int main() 
{
  // std::ofstream ofs(biasfilename);
  // std::ofstream ofs("biasdata/biasmap" + std::string(".bin"), std::ios::out | std::ios::binary);
  // if (ofs.fail()) {
  //   std::cout << "The bias file couldn't be opened. 'mkdir biasdata'" << std::endl;
  //   return -1;
  // }

  std::vector<std::unique_ptr<std::ofstream>> ofs_vector;
  for (int i = 0; i < totalnoiseNo; ++i) {
    std::string filename = biasfilename + "_" + std::to_string(i) + ".dat";
    ofs_vector.emplace_back(std::make_unique<std::ofstream>(filename));
        
    if (!ofs_vector.back()->is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return 1;
    }
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
  int divstep = int(totalstep/totalnoiseNo);
  int modstep = int(totalstep%totalnoiseNo);
  for (int l=0; l<totalnoiseNo+1; l++) {
    std::vector<std::vector<double>> biasdata;
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
        biasdata[i] = biaslist_fftw((i+l*divstep)*dN); // biaslist((i+l*divstep)*dN); // 
      }
// #ifdef _OPENMP
// #pragma omp critical
// #endif
      {
        count++;
	std::cout << "\rBiasGenerating : " << std::setw(3) << 100*count/divstep << "%" << std::flush;
      }
    }
    std::cout << std::endl;
    
    // for (size_t n=0; n<biasdata.size(); n++) {

      // for (size_t i=0; i<biasdata[0].size(); i++) {
      //   if (l<totalnoiseNo || n<modstep) ofs << biasdata[n][i] << ' ';
      // }
      // if (l<totalnoiseNo || n<modstep) ofs << std::endl;

      // ofs.write(reinterpret_cast<const char*>(biasdata[n].data()), biasdata[n].size()*sizeof(double));

      // std::cout << "\rExporting :      " << std::setw(3) << 100*n/biasdata.size() << "%" << std::flush;
    // }

    for (size_t n=0; n<biasdata.size(); n++) {
      int divcount=0;
      for (size_t m=0; m<biasdata[0].size(); m++) {
        *ofs_vector[divcount] << biasdata[n][m] << ' ';
        if ((m+1)%NL==0) {
          *ofs_vector[divcount] << std::endl;
          divcount++;
        }
      }
      std::cout << "\rExporting : " << std::setw(3) << 100*n/biasdata[0].size() << "%" << std::flush;
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

  {
    std::vector<std::vector<std::vector<std::complex<double>>>> biaslattice = fft_fftw(bk);

    LOOP{
      biaslist[i*NLnoise*NLnoise + j*NLnoise + k] = biaslattice[i][j][k].real();
    }
  }

  return biaslist;
}
