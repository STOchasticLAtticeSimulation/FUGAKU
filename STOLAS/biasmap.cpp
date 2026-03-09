#include "src/bias_noise.hpp"

std::vector<double> biaslist(double N);

int main() 
{
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------------

  for (int n = 0; n < internumber+2; n++) {
    int biasstepnum = firststep;
    if (n>0) biasstepnum = itpstep;
    if (n==internumber+1) biasstepnum = totalstep-firststep;

    std::vector<std::unique_ptr<std::ofstream>> ofs_vector;
    for (int i = 0; i < totalnoiseNo; i++) {
      std::string filename = biasfiledir + noisefilenamediv + std::to_string(i) + "_" + std::to_string(n) + ".bin";
      if (n==internumber+1) filename = biasfiledir + noisefilenamediv + std::to_string(i) + "_" + std::to_string(999) + ".bin";

      ofs_vector.emplace_back(std::make_unique<std::ofstream>(filename, std::ios::out | std::ios::binary));
          
      if (!ofs_vector.back()->is_open()) {
          std::cerr << "Failed to open file: " << filename << std::endl;
          return 1;
      }
    }

fftw_init_threads();
#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    std::cout << "Box size : " << NLnoise << "  Step number : " << biasstepnum << std::endl;
    
    int divstep = int((biasstepnum)/totalnoiseNo);
    int modstep = int((biasstepnum)%totalnoiseNo);
    
    for (int l=0; l<totalnoiseNo+1; l++) {
      std::vector<std::vector<double>> biasdata;
      if (l<totalnoiseNo) {
        biasdata = std::vector<std::vector<double>>(divstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
      } else {
        biasdata = std::vector<std::vector<double>>(modstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
      }

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i=0; i<divstep; i++) {
        int Nstep = i+l*divstep;
        if (n>0) Nstep += firststep-itpstep;
        if (n==internumber+1) Nstep += itpstep;
        if (l<totalnoiseNo || i<modstep) biasdata[i] = biaslist(Nstep*dN);
      }

      for (size_t n=0; n<biasdata.size(); n++) {
        int divcount=0;
        for (size_t m=0; m<biasdata[0].size(); m++) {
          ofs_vector[divcount]->write(reinterpret_cast<const char*>(&biasdata[n][m]), sizeof(double));
          if ((m+1)%NL==0) {
            divcount++;
          }
        }
      }

      std::cout << "\rBiasGenerating : " << std::setw(3) << 100*l/(totalnoiseNo) << "%" << std::flush;
    }
    std::cout << "\rBiasGenerating : 100%" << std::endl;
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

  std::vector<std::vector<std::vector<std::complex<double>>>> biaslattice = fft_fftw(bk);
  LOOP{
    biaslist[i*NLnoise*NLnoise + j*NLnoise + k] = biaslattice[i][j][k].real();
  }

  return biaslist;
}
