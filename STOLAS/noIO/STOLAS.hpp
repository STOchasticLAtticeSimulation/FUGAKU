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
#include <random>
#include <sys/time.h>
#include <complex>
#include <bit>
#include <boost/numeric/odeint.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

constexpr double euler_gamma = 0.57721566490153286061;
constexpr double LOG2 = 0.69314718055994530942;
const std::complex<double> II(0, 1);
std::normal_distribution<> dist(0., 1.);

constexpr int log2_int(unsigned int x) {
  return std::countr_zero(x);
}


#include "parameters.hpp"
#include "model.hpp"
#include "src/vec_op.hpp"
#include "src/array_op.hpp"
#include "src/fft.hpp"
#include "src/util.hpp"

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
std::array<double,NLnoiseAll> Ndata{};

// state_type phii{};
std::array<state_type,NLnoiseAll> Phidata{};
std::array<state_type,NLnoiseAll> phievol{};
std::array<state_type,NLnoiseAll> PhidataAv{};
// std::vector<std::vector<std::vector<double>>> phianimation;

std::array<double,NLnoiseAll> Nnoise{};
std::array<double,NLnoiseAll> Ntotal{};

std::array<std::array<double,NLnoiseAll>,NFIELDS+1> biaslist{};
std::array<std::array<double,NLnoiseAll>,NFIELDS+1> dwlist{};

// output
std::array<double,imax> disc_power{};
std::array<double,100*NLnoise> weightlist{};
std::array<std::array<double,NLnoise/2>,2> zetar{};
std::array<double,NLnoise/2> dzetar{};

#if MODEL==1
std::array<double,NLnoiseAll> N0list{};
std::array<bool,NLnoiseAll> brokenlist{};
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

  LOOP{
    phievol[NLnoise*NLnoise*i + NLnoise*j + k] = phii;
  }
}


// judge if point is in nsigma sphere shell
bool innsigma(int nx, int ny, int nz, int Num, double nsigma, double dn) {
  int nxt = (nx<=Num/2 ? nx : nx-Num);
  int nyt = (ny<=Num/2 ? ny : ny-Num);
  int nzt = (nz<=Num/2 ? nz : nz-Num);

  return std::abs(sqrt(nxt*nxt + nyt*nyt + nzt*nzt) - nsigma) <= dn/2.;
}

// judge real point
bool realpoint(int nx, int ny, int nz, int Num) {
  return (nx==0||nx==Num/2) && (ny==0||ny==Num/2) && (nz==0||nz==Num/2);
}

// judge independent complex point
bool complexpoint(int nx, int ny, int nz, int Num) {
  int nxt = (nx<=Num/2 ? nx : nx-Num);
  int nyt = (ny<=Num/2 ? ny : ny-Num);
  int nzt = (nz<=Num/2 ? nz : nz-Num);

  return (1<=nxt && nxt!=Num/2 && nyt!=Num/2 && nzt!=Num/2) ||
    (nxt==Num/2 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) || (nxt!=Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt!=Num/2) ||
    (nxt==0 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) ||
    (nxt==Num/2 && nyt==Num/2 && 1<=nzt && nzt!=Num/2) || (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt==Num/2) ||
    (nxt==0 && 1<=nyt && nyt!=Num/2 && nzt==0) ||
    (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==0) || (1<=nxt && nxt!=Num/2 && nyt==0 && nzt==Num/2) || (nxt==0 && nyt==Num/2 && 1<=nzt && nzt!=Num/2);
}

void dwlist_gen(double N, std::mt19937& engine) {
  int count = 0;
  double nsigma = sigma*exp(N);

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
    for (int i = 0; i < NLnoiseAll; i++) {
      dwlist[0][i] = out[i][0];
    }
    fftw_free(in);
    fftw_free(out);
    return;
  }

  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] /= sqrt(count);
    in[i][1] /= sqrt(count);
  }

  fftw_plan plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);

  for (int i = 0; i < NLnoiseAll; i++) {
    dwlist[0][i] = out[i][0];
  }

  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);
}

void biaslist1D(double N) {
  int count = 0;
  double nsigma = sigma*exp(N);
  
  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] = 0.0;
    in[i][1] = 0.0;
  }

  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      int idx = i*NLnoise*NLnoise + j*NLnoise + k;
      in[idx][0] = 1.0;
      in[idx][1] = 0.0;
      count++;
    }
  }

  if (count==0) {
    for (int i = 0; i < NLnoiseAll; i++) {
      biaslist[0][i] = out[i][0];
    }
    fftw_free(in);
    fftw_free(out);
    return;
  }
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] /= count;
  }

  fftw_plan plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);

  for (int i = 0; i < NLnoiseAll; i++) {
    biaslist[0][i] = out[i][0];
  }

  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);
}


void evolution(int seed) {
  double N = 0;
  int animationcount = 0; // for animation
  std::mt19937 engine(seed);

  for (size_t n=0; n<totalstep; n++){
    dwlist_gen(n*dN,engine);
    biaslist1D(n*dN);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<NLnoiseAll; i++){
      boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;
      state_type phi = phievol[i];
      
      #if MODEL==1
        double phiamp = sqrt(calPphi(N,phi,N0list[i],brokenlist[i]));
        double piamp = sqrt(calPpi(N,phi,N0list[i],brokenlist[i]));
        double crosscor = RecalPphipi(N,phi,N0list[i],brokenlist[i]);
      #elif MODEL==2
        double phiamp = sqrt(calPphi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]));
        double piamp = sqrt(calPpi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]));
        double crosscor = RecalPphipi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
      #else
        double phiamp = sqrt(calPphi(phi));
      #endif

      double dw = dwlist[0][i];
      if(i==0 && sweight) weightlist[n] = dw; // save weight data

      double Bias = biaslist[0][i];
      double GaussianFactor = 1./dNbias/sqrt(2*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);

      #if MODEL==2
        double Nstep = N;
        for (int dn=0;dn<(int)divdN;dn++) {
          stepper_noise.do_step(dphidN, phi, Nstep, dN/divdN);
          Nstep += dN/divdN;
        }
      #else
        stepper_noise.do_step(dphidN, phi, N, dN);
      #endif

      phi[0] += phiamp * dw * sqrt_dN;
      phi[0] += phiamp * bias * Bias * GaussianFactor * dN;
      phievol[i] = phi;

      #if MODEL==1
        if (crosscor > 0) {
          phi[1] += piamp * dw * sqrt_dN;
          phi[1] += piamp * bias * Bias * GaussianFactor * dN;
        } else {
          phi[1] -= piamp * dw * sqrt_dN;
          phi[1] -= piamp * bias * Bias * GaussianFactor * dN;
        }

        if (!brokenlist[i] && phi[0] < 0) {
          brokenlist[i] = true;
          N0list[i] = N;
        }
      #elif MODEL==2
        if (!broken1list[i] && phi[0] < phi1) {
          broken1list[i] = true;
          N1list[i] = N;
          if(i==0) std::cout << std::endl << "N1 = " << N << std::endl;
        }
        if (!broken2list[i] && phi[0] < phi2) {
          broken2list[i] = true;
          N2list[i] = N;
          if(i==0) std::cout << std::endl << "N2 = " << N << std::endl;
        }
      #endif

      if(strajectory && i==0){
        double Nd = N;
        save_trajectory(phi, Nd);
      }

      // animationcount++;
      // if(sanimation && animationcount>aninum-1){
      //   int ef = int(n/aninum);
      //   NFLOOP{
      //     phianimation[ef][i][2*nf] = phi[2*nf];
      //   }
      //   animationcount=0;
      // }
    }

    N += dN;
    std::cout << "\rLatticeSimulation   : " << std::setw(1) << int(100.*n/(double)totalstep) << "%" << std::flush;
  }
}


void evolutionNoise(int NoiseNo) {

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NLnoiseAll; i++) {
    double N = 0;
    int LatticePoint = i + NLnoiseAll*NoiseNo;
    
    std::mt19937 engine(1 + LatticePoint); // generate random noise for each thread
    std::normal_distribution<> dist(0., 1.);

    state_type phi = phievol[LatticePoint];

    #if MODEL==1
      double N0 = N0list[LatticePoint];
      bool broken = (brokenlist[LatticePoint]==0 ? false : true);
    #elif MODEL==2
      double N1 = N1list[LatticePoint];
      double N2 = N2list[LatticePoint];
      bool broken1 = (broken1list[LatticePoint]==0 ? false : true);
      bool broken2 = (broken1list[LatticePoint]==0 ? false : true);
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

    phievol[LatticePoint] = phi;
    Nnoise[LatticePoint] = N;
  }
}


void dNmap(int InterpolatingNo) {
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NLnoiseAll; i++) {
    int numstep = 0;
    int LatticePoint = i;
    double N = dN*(totalstep + InterpolatingNo*itpstep) + Nnoise[LatticePoint];
    state_type phi = phievol[LatticePoint];

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

      if(strajectory && i==0 && superH) {
        save_trajectory(phi, N);
      }
    }

    Ndata[LatticePoint] = N;
    phievol[LatticePoint] = phi;
  }
}


bool checkNfilefail() {
  return Nfilefail;
}


#endif
