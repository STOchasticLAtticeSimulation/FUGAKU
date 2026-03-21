#include "src/bias_noise.hpp"

std::vector<double> biaslist1D(double N);

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

    int divstep = int((biasstepnum)/totalnoiseNo);
    int modstep = int((biasstepnum)%totalnoiseNo);
    std::vector<std::vector<double>> biasdatadiv(divstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
    std::vector<std::vector<double>> biasdatamod(modstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));

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

    for (int l=0; l<totalnoiseNo+1; l++) {
      if (l<totalnoiseNo) {
        for (int i=0; i<divstep; i++) {
          int Nstep = i+l*divstep;
          if (n>0) Nstep += firststep-itpstep;
          if (n==internumber+1) Nstep += itpstep;
          if (l<totalnoiseNo || i<modstep) biasdatadiv[i] = biaslist1D(Nstep*dN);
        }

        for (size_t n = 0; n < biasdatadiv.size(); n++) {
          int divcount = 0;
          for (size_t m = 0; m < biasdatadiv[n].size(); m += NL) {
            ofs_vector[divcount]->write(
              reinterpret_cast<const char*>(&biasdatadiv[n][m]),
              sizeof(double) * NL
            );
            divcount++;
          }
        }
      } else {
        for (int i=0; i<divstep; i++) {
          int Nstep = i+l*divstep;
          if (n>0) Nstep += firststep-itpstep;
          if (n==internumber+1) Nstep += itpstep;
          if (l<totalnoiseNo || i<modstep) biasdatamod[i] = biaslist1D(Nstep*dN);
        }

        for (size_t n = 0; n < biasdatamod.size(); n++) {
          int divcount = 0;
          for (size_t m = 0; m < biasdatamod[n].size(); m += NL) {
            ofs_vector[divcount]->write(
              reinterpret_cast<const char*>(&biasdatamod[n][m]),
              sizeof(double) * NL
            );
            divcount++;
          }
        }
      }

      // std::cout << "\rBiasGenerating : " << std::setw(3) << 100*l/(totalnoiseNo) << "%" << std::flush;
    }
    // std::cout << "\rBiasGenerating : 100%" << std::endl;
  }

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}

std::vector<double> biaslist1D(double N) {
  int count = 0;
  double nsigma = sigma*exp(N);
  std::vector<double> biaslist(NLnoise*NLnoise*NLnoise,0);

  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] = 0.0;
    in[i][1] = 0.0;
  }

  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      int idx = i*NLnoise*NLnoise + j*NLnoise + k;
      in[idx][0] = 1.0;
      in[idx][1] = 0.0;
      count++;
    }
  }

  if (count==0) {
    fftw_free(in);
    fftw_free(out);
    return biaslist;
  }
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] /= count;
  }

  fftw_plan plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);

  for (int i = 0; i < NLnoiseAll; i++) {
    biaslist[i] = out[i][0];
  }

  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);

  return biaslist;
}
