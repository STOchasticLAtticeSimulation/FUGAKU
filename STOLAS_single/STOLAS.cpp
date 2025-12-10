#include "STOLAS.hpp"


// main function
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

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif

  int noisefiledirNo = atoi(argv[1]);

  Nfile.open(Nfileprefix + std::to_string(NLnoise) + std::string("_") + std::to_string(noisefiledirNo) + std::string(".dat"));
  logwfile.open(logwfileprefix + std::to_string(NLnoise) + std::string(".dat"), std::ios::app);
  Nfilefail = Nfile.fail();
  Nfile << std::setprecision(10);
  
  Ndata = std::vector<double>(NLnoise*NLnoise*NLnoise,0);

  if (sanimation) {
      Hdata = std::vector<std::vector<double>>(totalstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
      pidata = std::vector<std::vector<double>>(totalstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
    }
  
  for (int noiseNo = 0; noiseNo < totalnoiseNo; noiseNo++)
  {
    STOLAS(noisefiledirNo,noiseNo);

    if (!checknoisefile()) {
      std::cout << "The noise file couldn't be opened." << std::endl;
      return -1;
    }

    if (!checkbiasfile()) {
      std::cout << "The bias file couldn't be opened." << std::endl;
      return -1;
    }

    if (!noisebiassize()) {
      std::cout << "The box sizes of the noise and the bias are inconsistent." << std::endl;
      return -1;
    }

    if (checkNfilefail()) {
      std::cout << "The export file couldn't be opened. 'mkdir data'" << std::endl;
      return -1;
    }

    if(szeta) dNmap(noiseNo);
    if(sweight || noiseNo==0) weight();
  }

  if(spower) spectrum(Ndata,noisefiledirNo);
  if(sanimation) animation(noisefiledirNo);
  if(scompaction) compaction(Ndata,noisefiledirNo);

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  std::cout << std::endl;
  // -------------------------------------
}
