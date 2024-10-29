#include "bias_noise.hpp"

std::vector<double> biaslist(double N);

int main() 
{
  std::ofstream ofs(biasfilename);
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

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif

  std::cout << "Box size : " << NL << std::endl;
  
  int totalstep = ceil(log((NL/2-1)/sigma)/dN), count = 0;
  // std::vector<std::vector<double>> biasdata(totalstep, std::vector<double>(NL*NL*NL,0));
  int divnumber = 10;
  int divstep = int(totalstep/divnumber);
  int modstep = int(totalstep%divnumber);
  for (int l=0; l<divnumber+1; l++) {
    std::vector<std::vector<double>> biasdata;
    if (l<divnumber) {
      biasdata = std::vector<std::vector<double>>(divstep, std::vector<double>(NL*NL*NL,0));
    } else {
      biasdata = std::vector<std::vector<double>>(modstep, std::vector<double>(NL*NL*NL,0));
    }
    int count = 0;
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<divstep; i++) {
      if (l<divnumber || i<modstep) biasdata[i] = biaslist(i*dN);
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        count++;
	std::cout << "\rBiasGenerating : " << std::setw(3) << 100*count/divstep << "%" << std::flush;
      }
    }
    std::cout << std::endl;
    
    for (size_t n=0; n<biasdata.size(); n++) {
      /*
      for (size_t i=0; i<biasdata[0].size(); i++) {
        if (l<divnumber || n<modstep) ofs << biasdata[n][i] << ' ';
      }
      if (l<divnumber || n<modstep) ofs << std::endl;
      */
      ofs.write(reinterpret_cast<const char*>(biasdata[n].data()), biasdata[n].size()*sizeof(double));
      std::cout << "\rExporting : " << std::setw(3) << 100*n/biasdata.size() << "%" << std::flush;
    }
    std::cout << "\rExporting : 100%" << std::endl;
  }

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}


std::vector<double> biaslist(double N) {
  std::vector<std::vector<std::vector<std::complex<double>>>> bk(NL, std::vector<std::vector<std::complex<double>>>(NL, std::vector<std::complex<double>>(NL, 0)));
  int count = 0;
  double nsigma = sigma*exp(N);
  std::vector<double> biaslist(NL*NL*NL,0);
  
  LOOP{
    if (innsigma(i,j,k,NL,nsigma,dn)) {
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
    biaslist[i*NL*NL + j*NL + k] = biaslattice[i][j][k].real();
  }

  return biaslist;
}
