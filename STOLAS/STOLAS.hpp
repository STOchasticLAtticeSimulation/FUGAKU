#ifndef INCLUDED_STOLAS_
#define INCLUDED_STOLAS_

#define _USR_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <complex>
#include <sys/time.h>

#include <boost/numeric/odeint.hpp>

#define euler_gamma 0.57721566490153286061

#include "parameters.hpp"
#include "model.hpp"
#include "src/vec_op.hpp"
#include "src/array_op.hpp"
#include "src/fft.hpp"
#include "src/util.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


#include <random>
std::mt19937 engine(1); // Fixed seed
// std::random_device seed; // random seed
// std::mt19937 engine(seed());
std::normal_distribution<> dist(0., 1.);


const std::string Nfileprefix = sdatadir + "/" + model + "/Nmap_";
const std::string fieldfileprefix = sdatadir + "/" + model + "/field_";
const std::string trajectoryfileprefix = sdatadir + "/" + model + "/trajectory_";
const std::string animationfileprefix = sdatadir + "/" + model + "/animation/animation_";
const std::string powfileprefix = sdatadir + "/" + model + "/power_";
const std::string powsfileprefix = sdatadir + "/" + model + "/powers";
const std::string cmpfileprefix = sdatadir + "/" + model + "/compaction_";
const std::string prbfileprefix = sdatadir + "/" + model + "/probabilities";
const std::string logwfileprefix = sdatadir + "/" + model + "/logw_";
bool noisefilefail, biasfilefail, Nfilefail, sanisuperH = false;

int Nn, noisefiledirNo, noisefileNo;
std::ofstream Nfile, fieldfile, fieldfileA, trajectoryfile, powfile, powsfile, cmpfile, prbfile, logwfile;
std::vector<std::vector<double>> noisedataTr;
std::vector<std::vector<double>> biasdataTr;
std::vector<std::vector<std::vector<double>>> noisedata(NFIELDS+1);
std::vector<std::vector<std::vector<double>>> biasdata(NFIELDS+1);
std::array<double,NLnoiseAll> Ndata{};
std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D;

// state_type phii{};
std::array<state_type,NLnoiseAll> Phidata{};
std::array<state_type,NLnoiseAll> phievol{};
std::array<state_type,NLnoiseAll> PhidataAv{};
std::vector<std::vector<std::vector<double>>> phianimation;

std::array<double,NLnoiseAll> Nnoise{};
std::array<double,NLnoiseAll> Naverage{};
std::array<double,NLnoiseAll> Ntotal{};

#include "src/output.hpp"
#include "src/zoom.hpp"
#include "src/laplacian.hpp"


// -- functions -----------------------
void initialize(){
  phii[0] = PHI_INIT;
  phii[1] = DPHI_INIT;

  #if MODEL==3
  NFLOOP{
    phii[2*nf] = PSI_INIT;
    phii[2*nf+1] = DPSI_INIT;
  }
  #endif

  LOOPLONG{
    phievol[NLnoise*NLnoise*i + NLnoise*j + k] = phii;
  }
}


void STOLAS(int NoisefileDirNo, int NoisefileNo, int InterpolatingNo) {
  std::ifstream noisefile, biasfile;
  
  #if MODEL==3
  NFLOOP{
    noisefile.open(noisefiledir + std::to_string(NoisefileDirNo+nf-1) + noisefilenamediv + std::to_string(NoisefileNo) + "_" + std::to_string(InterpolatingNo) + std::string(".bin"), std::ios::binary);
    noisefilefail = noisefile.fail();
    
    if (!noisefile.fail()) {
      while (true) {
        std::vector<double> vv(NL);
        noisefile.read(reinterpret_cast<char*>(vv.data()), NL * sizeof(double));

        if (noisefile.gcount() != NL * sizeof(double)) {
          break;
        }
        noisedataTr.push_back(vv);
      }
      noisedata[nf-1] = transpose(noisedataTr);
      noisedataTr.clear();
      noisefile.close();
    }
  }
  #else
    noisefile.open(noisefiledir + std::to_string(NoisefileDirNo) + noisefilenamediv + std::to_string(NoisefileNo) + "_" + std::to_string(InterpolatingNo) + std::string(".bin"), std::ios::binary);
    noisefilefail = noisefile.fail();
    biasfile.open(biasfiledir + noisefilenamediv + std::to_string(NoisefileNo) + "_" + std::to_string(InterpolatingNo) + std::string(".bin"), std::ios::binary);
    biasfilefail = biasfile.fail();
      
    if (!noisefile.fail()) {
      while (true) {
        std::vector<double> vv(NL);
        noisefile.read(reinterpret_cast<char*>(vv.data()), NL * sizeof(double));

        if (noisefile.gcount() != NL * sizeof(double)) {
          break;
        }
        noisedataTr.push_back(vv);
      }
      noisedata[0] = transpose(noisedataTr);
      noisedataTr.clear();
      noisefile.close();
    }

    if (!biasfile.fail()) {
      while (true) {
        std::vector<double> vv(NL);
        biasfile.read(reinterpret_cast<char*>(vv.data()), NL * sizeof(double));

        if (biasfile.gcount() != NL * sizeof(double)) {
          break;
        }
        biasdataTr.push_back(vv);
      }
      biasdata[0] = transpose(biasdataTr);
      biasdataTr.clear();
      biasfile.close();
    }
  #endif

}


void evolution(int NoiseNo) {

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL; i++) {
    double N = 0;
    int numstep = 0;
    int animationcount = 0; // for animation
    #if MODEL==1 || MODEL==2
      double N0;
      bool broken = false;
    #endif
    
    state_type phi = phievol[i + NL*NoiseNo];

    // stepper
    boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;

    for (size_t n=0; n<noisedata[0][0].size(); n++) {
      #if MODEL==1 || MODEL==2
        double phiamp = sqrt(calPphi(N,phi,N0,broken));
        double piamp = sqrt(calPpi(N,phi,N0,broken));
        double crosscor = RecalPphipi(N,phi,N0,broken);
      #else
        double phiamp = sqrt(calPphi(phi));
      #endif

      #if MODEL==2
        for (int dn=0;dn<(int)divdN;dn++) {
          stepper_noise.do_step(dphidN, phi, N, dN/divdN);
          N += dN/divdN;
        }
      #else
        stepper_noise.do_step(dphidN, phi, N, dN);
        N += dN;
      #endif
  
      double dw = 0.;
      #if MODEL==3
        NFLOOP{
          dw = noisedata[nf-1][i][n];// dist(engine);// 
          phi[2*nf] += phiamp * dw * sqrt_dN;
        }
      #else
        dw = noisedata[0][i][n];
        double Bias = biasdata[0][i][n];
        double GaussianFactor = 1./dNbias/sqrt(2*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
        
        phi[0] += phiamp * dw * sqrt_dN;
        phi[0] += phiamp * bias * Bias * GaussianFactor * dN;
      #endif

      if(strajectory && i==0 && NoiseNo==0){
        double Nd = N;
        if(sanisuperH) Nd += dN*firststep;
        save_trajectory(phi, Nd);
      }

      animationcount++;
      if(sanimation && animationcount>aninum-1){
        int ef = int(n/aninum);
        if(sanisuperH) ef += int(firststep/aninum);
        NFLOOP{
          phianimation[ef][i + NL*NoiseNo][2*nf] = phi[2*nf];
        }
        animationcount=0;
      }
    }

    phievol[i + NL*NoiseNo] = phi;
  }
}


void evolutionNoise(int NoiseNo) {

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL; i++) {
    double N = 0;
    state_type phi = phievol[i + NL*NoiseNo];

    // stepper
    // boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;
    
    // while (EoN(phi)>0) {
    //   #if MODEL==3
    //     double psiamp = sqrt(calPpsi(phi));
    //     stepper_noise.do_step(dphidN, phi, N, dN);
    //     N += dN;
    //     NFLOOP{
    //       double dw = dist(engine);
    //       phi[2*nf] += psiamp * dw * sqrt_dN;
    //     }

    //     if(strajectory && i==0 && NoiseNo==0 && sanisuperH) {
    //       save_trajectory(phi, N+dN*totalstep);
    //     }
    //   #else
    //     double phiamp = sqrt(calPphi(phi));
    //     stepper_noise.do_step(dphidN, phi, N, dN);
    //     N += dN;
        
    //     double dw = dist(engine);
    //     phi[0] += phiamp * dw * sqrt_dN;

    //     if(strajectory && i==0 && NoiseNo==0 && sanisuperH) {
    //       save_trajectory(phi, N+dN*totalstep);
    //     }
    //   #endif
    // }

    // phievol[i + NL*NoiseNo] = phi;
    // Nnoise[i + NL*NoiseNo] = N;
  }
}


void dNmap(int NoiseNo, int InterpolatingNo) {
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL; i++) {
    int numstep = 0;
    double N = dN*(totalstep + InterpolatingNo*itpstep) + Nnoise[i + NL*NoiseNo];
    state_type phi = phievol[i + NL*NoiseNo];

    // Find zero crossing time
    typedef boost::numeric::odeint::runge_kutta_dopri5<state_type>base_stepper_type;
    auto stepper = make_dense_output(1.0e-15, 1.0e-15, base_stepper_type());

    stepper.initialize(phi, N, dN);
    state_type phil{};
    state_type phir{};
    state_type phim{};

    while (true){
      stepper.do_step(dphidN);
      phi = stepper.current_state();
      N = stepper.current_time();
      
      if (EoI(phi)<0){
        double precphi = 1.e+2;
        double Nl = N;
        double Nr = N - stepper.current_time_step();
        double Nmid = 0;

        while (precphi>Nprec){
          stepper.calc_state(Nl, phil);
          stepper.calc_state(Nr, phir);
          Nmid = (Nl*EoI(phir) - Nr*EoI(phil)) / (EoI(phir)-EoI(phil));
          stepper.calc_state(Nmid, phim);
          if (EoI(phim)>0){
            Nr = Nmid;
          }
          else {
            Nl = Nmid;
          }
          precphi = std::fabs(EoI(phim));
        }

        N = Nmid;
        break;
      }

      if(strajectory && i==0 && NoiseNo==0 && sanisuperH) {
        save_trajectory(phi, N);
      }
    }

    Ndata[i + NL*NoiseNo] = N;
    phievol[i + NL*NoiseNo] = phi;
  }
}


bool checknoisefile() {
  return noisefilefail;
}

bool checkbiasfile() {
  return biasfilefail;
}

bool checkNfilefail() {
  return Nfilefail;
}


#endif
