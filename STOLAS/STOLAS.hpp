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
const std::complex<double> II(0, 1);

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
bool noisefilefail, biasfilefail, Nfilefail, superH = false;

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

#if MODEL==1
std::array<double,NLnoiseAll> N0list{};
std::array<double,NLnoiseAll> brokenlist{};
#elif MODEL==2
std::array<double,NLnoiseAll> N1list{};
std::array<double,NLnoiseAll> N2list{};
std::array<double,NLnoiseAll> broken1list{};
std::array<double,NLnoiseAll> broken2list{};
#endif

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
    double N = (superH ? firststep*dN : 0);
    // if(superH) N = firststep*dN;
    int numstep = 0;
    int animationcount = 0; // for animation
    #if MODEL==1
      double N0 = N0list[i + NL*NoiseNo];
      bool broken = (brokenlist[i + NL*NoiseNo]==0 ? false : true);
    #elif MODEL==2
      double N1 = N1list[i + NL*NoiseNo];
      double N2 = N2list[i + NL*NoiseNo];
      bool broken1 = (broken1list[i + NL*NoiseNo]==0 ? false : true);
      bool broken2 = (broken1list[i + NL*NoiseNo]==0 ? false : true);
    #endif
    
    state_type phi = phievol[i + NL*NoiseNo];

    // stepper
    boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;

    for (size_t n=0; n<noisedata[0][0].size(); n++) {
      #if MODEL==1
        double phiamp = sqrt(calPphi(N,phi,N0,broken));
        double piamp = sqrt(calPpi(N,phi,N0,broken));
        double crosscor = RecalPphipi(N,phi,N0,broken);
      #elif MODEL==2
        double phiamp = sqrt(calPphi(N,phi,N1,N2,broken1,broken2));
        double piamp = sqrt(calPpi(N,phi,N1,N2,broken1,broken2));
        double crosscor = RecalPphipi(N,phi,N1,N2,broken1,broken2);
      #else
        double phiamp = sqrt(calPphi(phi));
      #endif
  
      double dw = 0.;
      #if MODEL==3
        NFLOOP{
          dw = noisedata[nf-1][i][n];
          phi[2*nf] += phiamp * dw * sqrt_dN;
        }
      #else
        dw = noisedata[0][i][n];
        double Bias = biasdata[0][i][n];
        double GaussianFactor = 1./dNbias/sqrt(2*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);

        #if MODEL==2
          for (int dn=0;dn<(int)divdN;dn++) {
            stepper_noise.do_step(dphidN, phi, N, dN/divdN);
            N += dN/divdN;
          }
        #else
          stepper_noise.do_step(dphidN, phi, N, dN);
          N += dN;
        #endif
        
        phi[0] += phiamp * dw * sqrt_dN;
        phi[0] += phiamp * bias * Bias * GaussianFactor * dN;
      #endif

      #if MODEL==1
        if (crosscor > 0) {
          phi[1] += piamp * dw * sqrt_dN;
          phi[1] += piamp * bias * Bias * GaussianFactor * dN;
        } else {
          phi[1] -= piamp * dw * sqrt_dN;
          phi[1] -= piamp * bias * Bias * GaussianFactor * dN;
        }

        if (brokenlist[i + NL*NoiseNo]==0 && phi[0] < 0) {
          N0 = N;
          broken = true;
          brokenlist[i + NL*NoiseNo] = 1;
          N0list[i + NL*NoiseNo] = N;
        }
      #elif MODEL==2
        if (broken2list[i + NL*NoiseNo]==0 && phi[0] < phi1) {
          N1 = N;
          broken1 = true;
          broken1list[i + NL*NoiseNo] = 1;
          N1list[i + NL*NoiseNo] = N;
        }
        if (broken2list[i + NL*NoiseNo]==0 && phi[0] < phi2) {
          N2 = N;
          broken2 = true;
          broken2list[i + NL*NoiseNo] = 1;
          N2list[i + NL*NoiseNo] = N;
        }
      #endif

      if(strajectory && i==0 && NoiseNo==0){
        double Nd = N;
        if(superH) Nd += dN*firststep;
        save_trajectory(phi, Nd);
      }

      animationcount++;
      if(sanimation && animationcount>aninum-1){
        int ef = int(n/aninum);
        if(superH) ef += int(firststep/aninum);
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

    #if MODEL==1
      double N0 = N0list[i + NL*NoiseNo];
      bool broken = (brokenlist[i + NL*NoiseNo]==0 ? false : true);
    #elif MODEL==2
      double N1 = N1list[i + NL*NoiseNo];
      double N2 = N2list[i + NL*NoiseNo];
      bool broken1 = (broken1list[i + NL*NoiseNo]==0 ? false : true);
      bool broken2 = (broken1list[i + NL*NoiseNo]==0 ? false : true);
    #endif

    // stepper
    boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;
    
    while (EoN(phi)>0) {
      #if MODEL==3
        double psiamp = sqrt(calPpsi(phi));
        stepper_noise.do_step(dphidN, phi, N, dN);
        N += dN;
        NFLOOP{
          double dw = dist(engine);
          phi[2*nf] += psiamp * dw * sqrt_dN;
        }
      #elif MODEL==2
        double phiamp = sqrt(calPphi(N,phi,N1,N2,broken1,broken2));
        double piamp = sqrt(calPpi(N,phi,N1,N2,broken1,broken2));
        double crosscor = RecalPphipi(N,phi,N1,N2,broken1,broken2);
        
        stepper_noise.do_step(dphidN, phi, N, dN);
        N += dN;
        
        double dw = dist(engine);
        phi[0] += phiamp * dw * sqrt_dN;
      #elif MODEL==1
        double phiamp = sqrt(calPphi(N,phi,N0,broken));
        double piamp = sqrt(calPpi(N,phi,N0,broken));
        double crosscor = RecalPphipi(N,phi,N0,broken);
      #else
        double phiamp = sqrt(calPphi(phi));
        
        stepper_noise.do_step(dphidN, phi, N, dN);
        N += dN;
        
        double dw = dist(engine);
        phi[0] += phiamp * dw * sqrt_dN;
      #endif
      if(strajectory && i==0 && NoiseNo==0 && superH) {
        save_trajectory(phi, N+dN*totalstep);
      }
    }

    phievol[i + NL*NoiseNo] = phi;
    Nnoise[i + NL*NoiseNo] = N;
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

      if(strajectory && i==0 && NoiseNo==0 && superH) {
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
