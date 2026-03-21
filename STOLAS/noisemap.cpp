#include "src/bias_noise.hpp"

std::vector<double> dwlist(double N);

int main(int argc, char* argv[])
{
  if (argc!=2) {
    std::cout << "Specify the noise file number correctly." << std::endl;
    return -1;
  }
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------------

  for (int n = 0; n < internumber+2; n++) {
    int noisestepnum = firststep;
    if (n>0) noisestepnum = itpstep;
    if (n==internumber+1) noisestepnum = totalstep-firststep;

    int divstep = int((noisestepnum)/totalnoiseNo);
    int modstep = int((noisestepnum)%totalnoiseNo);
    std::vector<std::vector<double>> noisedatadiv(divstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
    std::vector<std::vector<double>> noisedatamod(modstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));

    std::vector<std::unique_ptr<std::ofstream>> ofs_vector;
    for (int i = 0; i < totalnoiseNo; i++) {
      std::string filename = noisefiledir + std::string(argv[1]) + noisefilenamediv + std::to_string(i) + "_" + std::to_string(n) + ".bin";
      if (n==internumber+1) filename = noisefiledir + std::string(argv[1]) + noisefilenamediv + std::to_string(i) + "_" + std::to_string(999) + ".bin";

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

    std::cout << "Box size : " << NLnoise << "  Step number : " << noisestepnum << std::endl;

    for (int l=0; l<totalnoiseNo+1; l++) {
      if (l<totalnoiseNo) {
        for (int i=0; i<divstep; i++) {
          int Nstep = i+l*divstep;
          if (n>0) Nstep += firststep-itpstep;
          if (n==internumber+1) Nstep += itpstep;
          if (l<totalnoiseNo || i<modstep) noisedatadiv[i] = dwlist(Nstep*dN);
        }

        for (size_t n = 0; n < noisedatadiv.size(); n++) {
          int divcount = 0;
          for (size_t m = 0; m < noisedatadiv[n].size(); m += NL) {
            ofs_vector[divcount]->write(
              reinterpret_cast<const char*>(&noisedatadiv[n][m]),
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
          if (l<totalnoiseNo || i<modstep) noisedatamod[i] = dwlist(Nstep*dN);
        }

        for (size_t n = 0; n < noisedatamod.size(); n++) {
          int divcount = 0;
          for (size_t m = 0; m < noisedatamod[n].size(); m += NL) {
            ofs_vector[divcount]->write(
              reinterpret_cast<const char*>(&noisedatamod[n][m]),
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

std::vector<double> dwlist(double N) {
  int count = 0;
  double nsigma = sigma*exp(N);
  std::vector<double> dwlist(NLnoise*NLnoise*NLnoise,0);

  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] = 0.0;
    in[i][1] = 0.0;
  }

  LOOP{
    int idx = i*NLnoise*NLnoise + j*NLnoise + k;
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      if (realpoint(i,j,k,NLnoise)) {
	in[idx][0] = dist(engine);
	count++;
      } else if (complexpoint(i,j,k,NLnoise)) {
	in[idx][0] = dist(engine)/sqrt(2);
	in[idx][1] = dist(engine)/sqrt(2);
	count++;
      }
    }
  }

  // reflection
  int ip, jp, kp; // reflected index
  LOOP{
    int idx = i*NLnoise*NLnoise + j*NLnoise + k;
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      if (!(realpoint(i,j,k,NLnoise)||complexpoint(i,j,k,NLnoise))) {
        ip = (i==0 ? 0 : NLnoise-i);
        jp = (j==0 ? 0 : NLnoise-j);
        kp = (k==0 ? 0 : NLnoise-k);
	
	in[idx][0] = in[ip*NLnoise*NLnoise + jp*NLnoise + kp][0];
	in[idx][1] = -in[ip*NLnoise*NLnoise + jp*NLnoise + kp][1];
	count++;
      }
    }
  }

  if (count==0) {
    fftw_free(in);
    fftw_free(out);
    return dwlist;
  }

  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] /= sqrt(count);
    in[i][1] /= sqrt(count);
  }

  fftw_plan plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);

  for (int i = 0; i < NLnoiseAll; i++) {
    dwlist[i] = out[i][0];
  }

  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);

  return dwlist;
}
