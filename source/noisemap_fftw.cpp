#include "bias_noise.hpp"

std::vector<double> dwlist(double N);
std::vector<double> dwlist_fftw(double N);

int main(int argc, char* argv[]) 
{
  if (argc!=2) {
    std::cout << "Specify the noise file number correctly." << std::endl;
    return -1;
  }

  // std::ofstream ofs(noisefilename + std::string(argv[1]) + //std::string(".dat"));
	// 	    std::string(".bin"), std::ios::out | std::ios::binary);//, std::ios::app);


  // std::ofstream ofs(noisefilename + std::string(argv[1]) + std::string(".dat"));
  // if (ofs.fail()) {
  //   std::cout << "The noise file couldn't be opened. 'mkdir noisedata'" << std::endl;
  //   return -1;
  // }

  std::vector<std::unique_ptr<std::ofstream>> ofs_vector;
  for (int i = 0; i < totalnoiseNo; ++i) {
    std::string filename = noisefilename + std::string(argv[1]) + "_" + std::to_string(i) + ".bin";
    ofs_vector.emplace_back(std::make_unique<std::ofstream>(filename, std::ios::out | std::ios::binary));
        
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

    // for (size_t n=0; n<noisedata.size(); n++) {
    //   ofs.write(reinterpret_cast<const char*>(noisedata[n].data()), noisedata[n].size()*sizeof(double));
      
    //   std::cout << "\rExporting :       " << std::setw(3) << 100*n/noisedata.size() << "%" << std::flush;
    // }

    for (size_t n=0; n<noisedata.size(); n++) {
      int divcount=0;
      for (size_t m=0; m<noisedata[0].size(); m++) {
        // *ofs_vector[divcount] << noisedata[n][m] << ' ';
        ofs_vector[divcount]->write(reinterpret_cast<const char*>(&noisedata[n][m]), sizeof(double));
        if ((m+1)%NL==0) {
          // *ofs_vector[divcount] << std::endl;
          divcount++;
        }
      }
      std::cout << "\rExporting : " << std::setw(3) << 100*n/noisedata[0].size() << "%" << std::flush;
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
