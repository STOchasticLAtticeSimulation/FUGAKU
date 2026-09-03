#include "STOLAS.hpp"

std::ifstream Ninfile;

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

  int noisefiledirNo = atoi(argv[1]);
  std::string NLfilename = std::to_string(NLnoise) + std::string("_") + std::to_string(NFIELDS) + std::string("_") + std::to_string(noisefiledirNo);
  std::string InterFileName = NLfilename + std::string("_") + std::to_string(0);
  std::string SigmaFileName = std::string("_") + std::to_string(int(calPzeta0));
  
  std::cout << Nfileprefix + InterFileName + SigmaFileName << std::endl;
  Ninfile.open(Nfileprefix + InterFileName + SigmaFileName + std::string(".bin"));
  if (Ninfile.fail()) std::cout << "Nfile fail" << std::endl;

  Ninfile.read(reinterpret_cast<char*>(Ndata.data()), sizeof(double) * NLnoiseAll);

  // if(spower) spectrum(Ndata,noisefiledirNo);
  mu2(Ndata, noisefiledirNo);
  if(scompaction) compaction(Ndata,noisefiledirNo);


  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  std::cout << std::endl;
  // -------------------------------------
}
