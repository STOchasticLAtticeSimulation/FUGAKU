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

  int seed_val = atoi(argv[1]);
  std::mt19937 engine(seed_val);

  std::cout << "Noise seed No. " << seed_val << std::endl;
  std::cout << "model : " << model << std::endl;

  // initialize fields
  initialize();

  // Solve fields
  for (size_t interpolatingnumber = 0; interpolatingnumber < internumber+1; interpolatingnumber++) {
    OpenFiles(seed_val, interpolatingnumber);
    
    if (checkNfilefail()) {
      std::cout << "The export file couldn't be opened. 'mkdir data'" << std::endl;
      return -1;
    }

    if(interpolatingnumber==0) evolution(seed_val,engine,0,firststep);
    else  evolution(seed_val,engine,firststep-itpstep,firststep);
    
    Phidata = phievol; // save field values for zoom
    evolution(seed_val,engine,firststep,totalstep);

    if(EoI_noise && MeanNumber!=0) {
      PhidataAv = phievol;
      for (int av = 0; av < MeanNumber; av++) {
        evolutionNoise(seed_val,av); // Adding the noise until EoN()
        dNmap(interpolatingnumber);
        Ntotal += Ndata;
        Ndata.fill(0.0); // reset vector
        phievol = PhidataAv;
        std::cout << "\rAveraging           : " << std::setw(1) << int(100.*(av+1)/(double)MeanNumber) << "%" << std::flush;
      }
      std::cout << std::endl;
      for (int i=0; i<NLnoiseAll; i++) Ndata[i] = Ntotal[i]/(double)MeanNumber;
    }
    else if(EoI_noise && MeanNumber==0) {
      evolutionNoise(seed_val,0); // Adding the noise until EoN()
      dNmap(interpolatingnumber);
    }
    else {
      dNmap(interpolatingnumber);
    }

    save_zeta(); // save delta N map
    if(spower) spectrum(Ndata,interpolatingnumber);
    if(sfield) save_field();
    if(sweight) weight(seed_val);
    if(scompaction) compaction(Ndata,seed_val);
    #if MODEL==2
      save_N1N2(seed_val);
    #endif
    Nfile.close();
    fieldfile.close();


    // Interpolation
    if (interpolatingnumber!=internumber) {
      std::vector<int> shift{0,0,0}; // std::vector<int> shift = findNMaxBox(Ndata);

      // phievol = Phidata; // reset field values
      Ntotal.fill(0.0); // reset vector
      InterpolatingPhi(shift);

      std::cout << "Number of interpolation: " << interpolatingnumber+1 << std::endl;
    }
  }

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------
}
