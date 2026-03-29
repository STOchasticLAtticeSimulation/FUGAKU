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
  // ---------------------------------

fftw_init_threads();
#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif

  int noisefiledirNo = atoi(argv[1]);
  int seed_val = atoi(argv[1]);

  std::cout << "Noise seed No. : " << noisefiledirNo << std::endl;
  std::cout << "model : " << model << std::endl;

  // initialize fields
  initialize();

  // Solve fields
  for (size_t interpolatingnumber = 0; interpolatingnumber < internumber+1; interpolatingnumber++) {
    OpenFiles(noisefiledirNo, interpolatingnumber);
    
    if (checkNfilefail()) {
      std::cout << "The export file couldn't be opened. 'mkdir data'" << std::endl;
      return -1;
    }


    evolution(seed_val);
    dNmap(interpolatingnumber);

    save_zeta(); // svave delta N map
    if(sfield) save_field();
    if(spower) spectrum(Ndata,interpolatingnumber);
    if(sweight) weight();
    if(scompaction) compaction(Ndata,noisefiledirNo);
    Nfile.close();
    fieldfile.close();


    // // Interpolation
    // if (interpolatingnumber!=internumber) {
    //   std::vector<int> shift{0,0,0}; // std::vector<int> shift = findNMaxBox(Ndata);

    //   phievol = Phidata; // reset field values
    //   Ntotal.fill(0.0); // reset vector
    //   InterpolatingPhi(shift);

    //   std::cout << "Number of interpolation: " << interpolatingnumber+1 << std::endl;
    // }
  }

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------
}
