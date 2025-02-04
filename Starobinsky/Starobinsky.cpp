#include "../source/STOLAS.hpp"

//--------- User may change ---------
// Model parameters
const std::string model = "Starobinsky"; // Name of the model
const std::vector<double> phii{0.0193,-5.45e-7}; // Initial conditions {field,derivative}
//-----------------------------------


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

  // file initialize
  std::ofstream Ninitfile(sdatadir + "/Nmap_" + std::to_string(NLnoise) + std::string("_") + std::to_string(noisefiledirNo) + std::string(".dat"));
  Ninitfile.close();
  std::ofstream logwinitfile(sdatadir + "/logw_" + std::to_string(NLnoise) + std::string("_") + std::to_string(noisefiledirNo) + std::string(".dat"));
  logwinitfile.close();
  
  for (int noiseNo = 0; noiseNo < totalnoiseNo; noiseNo++)
  {
    STOLAS stolas(model,dN,sourcedir,noisefiledirNo,phii,bias,Nbias,dNbias,noiseNo);

    if (!stolas.checknoisefile()) {
      std::cout << "The noise file couldn't be opened." << std::endl;
      return -1;
    }

    if (!stolas.checkbiasfile()) {
      std::cout << "The bias file couldn't be opened." << std::endl;
      return -1;
    }

    if (!stolas.noisebiassize()) {
      std::cout << "The box sizes of the noise and the bias are inconsistent." << std::endl;
      return -1;
    }

    if (stolas.Nfilefail()) {
      std::cout << "The export file couldn't be opened. 'mkdir data'" << std::endl;
      return -1;
    }

    if(szeta) stolas.dNmap(noiseNo);
    if(noiseNo==0) stolas.weight();
  }


  // if(spower) stolas.spectrum();
  // if(sanimation) stolas.animation();
  // if(scompaction) stolas.compaction();

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  std::cout << std::endl;
  // -------------------------------------
}
