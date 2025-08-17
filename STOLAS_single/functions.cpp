#include "STOLAS.hpp"

std::vector<std::vector<double>> Nindata;
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
  std::cout << Nfileprefix + std::to_string(NLnoise) + "_" + std::to_string(noisefiledirNo) << std::endl;
  Ninfile.open(Nfileprefix + std::to_string(NLnoise) + "_" + std::to_string(noisefiledirNo) + std::string(".dat"));
  if (Ninfile.fail()) std::cout << "Nfile fail" << std::endl;

  double first, second;
  while (Ninfile >> first >> second) {
    Nindata.push_back({first, second});
  }
  Ninfile.close();

  std::sort(Nindata.begin(), Nindata.end(), compareByFirstColumn);
  Ndata = std::vector<double>(NLnoise*NLnoise*NLnoise,0);

  for (size_t ii = 0; ii < Nindata.size(); ii++) {
    Ndata[ii] = Nindata[ii][1];
  }
  std::cout << "Sort done" << std::endl;

  if(spower) spectrum(Ndata,noisefiledirNo);
  if(scompaction) compaction(Ndata,noisefiledirNo);


  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  std::cout << std::endl;
  // -------------------------------------
}
