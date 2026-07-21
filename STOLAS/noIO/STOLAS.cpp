#include "STOLAS.hpp"

// main function
int main(int argc, char* argv[])
{
  if (argc!=2) {
    std::cout << "Specify the noise file number correctly." << std::endl;
    return -1;
  }
  
  if(NLnoise!=pow(2,NLpower)) {
    std::cout << "NLnoise and NLpower is not the same." << std::endl;
    return -1;
  }
  
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // ---------------------------------

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif

  init_fftw_global();
  int seed_val = atoi(argv[1]);
  std::mt19937 engine(seed_val);//engine(0);//
  std::mt19937 engine_int(seed_val);

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

    if(interpolatingnumber==0) {
      evolution(seed_val,engine,0,firststep,interpolatingnumber);
    }
    else {
      evolution(seed_val,engine,firststep-itpstep,firststep,interpolatingnumber);
    }
    
    Phidata = phievol; // save field values for zoom
    evolution(seed_val,engine,firststep,totalstep,interpolatingnumber);

    if(EoI_noise && MeanNumber!=0 && interpolatingnumber==internumber) {
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
    else if(EoI_noise && MeanNumber==0 && interpolatingnumber==internumber) {
      evolutionNoise(seed_val,0); // Adding the noise until EoN()
      dNmap(interpolatingnumber);
    }
    else {
      std::cout << "start dNmap" << std::endl;
      dNmap(interpolatingnumber);
    }

    // omit eternal inflation
    double Naverage = 0;
    double Ntyp = 50.;
    for (int i = 0; i < NLnoiseAll; i++) {
      if (Ndata[i] < Ntyp) Naverage += Ndata[i];
    }
    Naverage /= NLnoiseAll;

    double Nvariance = 0;
    for (int i = 0; i < NLnoiseAll; i++) {
      if (Ndata[i] < Ntyp) Nvariance += pw2(Ndata[i] - Naverage);
    }
    Nvariance /= NLnoiseAll;
    // std::cout << Naverage << ' ' << sqrt(Nvariance) << std::endl;

    double Nmax = Naverage + 10.*sqrt(Nvariance);
    for (int i = 0; i < NLnoiseAll; i++) {
       if (Ndata[i] > Nmax) Ndata[i] = Nmax;
    }

// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
//     LOOP {
//       int idx = index(i, j, k);
//       if (Ndata[idx] > 20.){
//         int ip1 = INCREMENT(i);
//         int im1 = DECREMENT(i);
//         int jp1 = INCREMENT(j);
//         int jm1 = DECREMENT(j);
//         int kp1 = INCREMENT(k);
//         int km1 = DECREMENT(k);

//         int neighbors[6] = {
//           index(ip1, j, k), index(im1, j, k),
//           index(i, jp1, k), index(i, jm1, k),
//           index(i, j, kp1), index(i, j, km1)
//         };

//         double sum = 0.0;
//         int count = 0;

//         for (int n = 0; n < 6; n++) {
//           if (Ndata[neighbors[n]] <= 20.) {
//             sum += Ndata[neighbors[n]];
//             count++;
//           }
//         }

//         if (count > 0) {
//           Ndata[idx] = sum / (double)count;
//         }
//       }
//     }

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
      int shiftx = (int)(NLnoise*0.25*(1.0+std::erf(dist(engine_int)*inv_sqrt2)));// NLnoise/4; //
      int shifty = (int)(NLnoise*0.25*(1.0+std::erf(dist(engine_int)*inv_sqrt2)));
      int shiftz = (int)(NLnoise*0.25*(1.0+std::erf(dist(engine_int)*inv_sqrt2)));
      shift = {shiftx,shifty,shiftz};
      std::cout << shiftx << ' ' << shifty << ' ' << shiftz << std::endl;
      // findNMaxBox(Ndata);
      Ntotal.fill(0.0); // reset vector
      InterpolatingPhi(shift);

      for (int i=0; i<NLnoiseAll; i++){
        #if MODEL==1
          if(phievol[i][0]>0) brokenlist[i] = false;
        #elif MODEL==2
          if (phievol[i][0] > phi1) broken1list[i] = false;
          if (phievol[i][0] > phi2) broken2list[i] = false;
        #endif
      }

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
