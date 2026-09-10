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
// #include <bit>
#include <boost/numeric/odeint.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

constexpr double euler_gamma = 0.57721566490153286061;
constexpr double LOG2 = 0.69314718055994530942;
const std::complex<double> II(0, 1);
std::normal_distribution<> dist(0., 1.);

// constexpr int log2_int(unsigned int x) {
//   return std::countr_zero(x);
// }


#include "parameters.hpp"
const double sqrt_dN = std::sqrt(dN);
constexpr int NLnoiseAll = NLnoise*NLnoise*NLnoise;
constexpr int NLnoiseHalf = NLnoise/2 + 1; // halved (fastest) dimension for r2c/c2r FFTW transforms
constexpr int NLnoiseHalfAll = NLnoise*NLnoise*NLnoiseHalf;
const int totalstep = ceil(log((NLnoise/2-1)/sigma)/dN); // Total number of time step with noise
const int firststep = ceil(log(nsigmareset/sigma)/dN);
const int itpstep = ceil((log(nsigmareset/sigma)-log(nsigmareset/sigma/2.))/dN);
constexpr double dx = LL/NLnoise; // Spacing of each lattice
// constexpr int NLpower = log2_int(NLnoise);
constexpr double imax_double = LOG2*(NLpower-1) / dlogn;
constexpr int imax = int(imax_double) + (imax_double > int(imax_double));
const double inv_sqrt2 = 1./sqrt(2.);
const double sig2 = sigma*sigma;
const double sig3 = sigma*sigma*sigma;


#include "model.hpp"
#include "src/vec_op.hpp" // use in laplacian
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
const std::string mu2fileprefix = sdatadir + "/" + model + "/mu2_";
const std::string k3fileprefix = sdatadir + "/" + model + "/k3_";

bool Nfilefail, superH = false;
bool FFTwisdomFirst = false;

// int noisefiledirNo, noisefileNo;
std::ofstream Nfile, fieldfile, fieldfileA, trajectoryfile, powfile, powsfile, cmpfile, prbfile, logwfile, Noisefile, mu2file, k3file;
std::array<double,NLnoiseAll> Ndata{};
std::array<double,NLnoiseAll> Nnoise{}; // use for EoN noise
std::array<double,NLnoiseAll> Ntotal{}; // use for averaging

std::array<state_type,NLnoiseAll> phievol{};
std::array<state_type,NLnoiseAll> PhidataAv{}; // use for averaging
std::array<state_type,NLnoiseAll> Phidata{}; // use for zoom


std::array<std::array<double,NLnoiseAll>,NFIELDS+1> biaslist{};
std::array<std::array<double,NLnoiseAll>,NFIELDS+1+1> dwlist{};
std::array<double,NLnoiseAll> laplacian{};
std::array<double,NLnoiseAll> mutwo{};
std::array<double,NLnoiseAll> kthree{};

// output
std::array<double,imax> disc_power{};
std::array<double,100*NLnoise> weightlist{};
std::array<double,100*NLnoise> weightbool{};
std::array<std::array<double,NLnoise/2>,2> zetar{};
std::array<double,NLnoise/2> dzetar{};

#if MODEL==1
  std::array<double,NLnoiseAll> N0list{};
  std::array<bool,NLnoiseAll> brokenlist{};
#elif MODEL==2
  std::array<double,NLnoiseAll> N1list{};
  std::array<double,NLnoiseAll> N2list{};
  std::array<bool,NLnoiseAll> broken1list{};
  std::array<bool,NLnoiseAll> broken2list{};
  const double dNsub = dN / divdN;
  const int ndiv  = (int)divdN;
#endif

#include "src/noise_bias.hpp"
#include "src/zoom.hpp"
#include "src/laplacian.hpp"
#include "src/output.hpp"
#if MODEL==2
#include "src/simd_rk4.hpp"
#endif


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

  // Parallelized with the same schedule(static) as the main per-point loops
  // so first-touch NUMA placement matches how the data is actually accessed
  // afterwards (matters on CMG/NUMA machines like Fugaku's A64FX; harmless
  // on a single-socket Mac).
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  LOOP{
    phievol[NLnoise*NLnoise*i + NLnoise*j + k] = phii;
  }
}


void evolution(int seed, int starttime, int endtime, int InterpolatingNo) {
  double N = (starttime + InterpolatingNo*itpstep)*dN;
  int animationcount = (starttime==firststep ? firststep%aninum : 0); // for animation
  int animationstep = (starttime==firststep ? firststep/((double)aninum) : 0);

  for (size_t n=starttime; n<endtime; n++){
    biaslist1D(n*dN);

    #if MODEL==3
      NFLOOP dwlist_gen(n*dN,seed,n,nf-1);
    #else
      dwlist_gen(n*dN,seed,n,0); // for phi
      dwlist_gen(n*dN,seed,n,1); // for pi
    #endif

    if(snoisemap){
      if(n==totalstep-500) {
        Noisefile.open(sdatadir + "/" + model + "/noisedata/map_" + std::to_string(n) + std::string(".bin"), std::ios::binary);
        Noisefile << std::setprecision(10);
        Noisefile.write(reinterpret_cast<const char*>(&dwlist[0]), sizeof(double) * NLnoiseAll);
        Noisefile.close();
      }
    }

#if MODEL==2
    // SIMD-batched drift step (NEON/SVE, see src/simd_rk4.hpp): this replaces
    // the ndiv-substep RK4 loop below, which profiling showed is >99% of the
    // per-step lattice cost. phievol is read-only here; results land in
    // driftedPhievol and are picked up per-point below, after calPphi/calPpi
    // /RecalPphipi (which need the pre-drift state).
    simd_rk4_drift_batch(phievol, driftedPhievol, ndiv, dNsub);
#endif

    // static: matches simd_rk4_drift_batch's schedule(static) above (so the
    // thread that just computed driftedPhievol[i] is also the one reading it
    // here -- keeps the two passes CMG/NUMA-local on Fugaku) and, since A64FX
    // cores are homogeneous (unlike this Mac's P/E core split, where guided's
    // dynamic rebalancing actually helps), avoids paying guided's per-chunk
    // dispatch overhead for no load-balancing benefit.
    // GaussianFactor depends only on N (the bias-window envelope), not on
    // the lattice point i -- hoisted out of the point loop, where it used to
    // be recomputed (1 exp() call) for all NLnoiseAll points every step.
    double GaussianFactor = 1./dNbias/sqrt(2.*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i=0; i<NLnoiseAll; i++){
      state_type phi = phievol[i];

      #if MODEL==1
        double phiamp = sqrt(calPphi(N,phi,N0list[i],brokenlist[i]));
        double piamp = sqrt(calPpi(N,phi,N0list[i],brokenlist[i]));
        double crosscor = RecalPphipi(N,phi,N0list[i],brokenlist[i]);
      #elif MODEL==2
        double calPphival = calPphi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double calPpival = calPpi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double crosscor = RecalPphipi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double phiamp = sqrt(calPphival);
      #else
        double phiamp = sqrt(calPphi(phi));
      #endif

      double dw = dwlist[0][i];
      if(i==0 && sweight){
        weightlist[n] = dw; // save weight data
        weightbool[n] = true;
      }

      double Bias = biaslist[0][i];

      #if MODEL==2
        double phi_old = phi[0]; // for reflective boundary
        phi[0] = driftedPhievol[i][0];
        phi[1] = driftedPhievol[i][1];
      #else
        boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;
        stepper_noise.do_step(dphidN, phi, N, dN);
      #endif

      #if MODEL==2
        double dwpi = dwlist[1][i];

        double b2 = crosscor*crosscor;
        double combi1 = sqrt(4.*b2 + pw2(calPphival-calPpival));

        double sqrtlam1 = sqrt(0.5*(calPphival+calPpival + combi1));
        double eig2 = 0.5*(calPphival+calPpival - combi1);
        if(eig2 < 1e-14*(calPphival+calPpival)) eig2 = 0.;

        double sqrtlam2 = sqrt(eig2);

        double denomplus = sqrt(4.*b2 + pw2(calPphival-calPpival + combi1));
        double denomminus = sqrt(4.*b2 + pw2(calPphival-calPpival - combi1));

        double vplus1 = (calPphival-calPpival + combi1)/denomplus;
        double vplus2 = 2.*crosscor/denomplus;
        double vminus1 = (calPphival-calPpival - combi1)/denomminus;
        double vminus2 = 2.*crosscor/denomminus;

        double biaseddw = dw + bias * Bias * GaussianFactor * sqrt_dN;

        phi[0] += (sqrtlam1*vplus1*biaseddw + sqrtlam2*vminus1*dwpi) * sqrt_dN;
        phi[1] += (sqrtlam1*vplus2*biaseddw + sqrtlam2*vminus2*dwpi) * sqrt_dN;

        // reflective boundary
        if (phi_old <= phi1 && phi[0] > phi1) {
          phi[0] = 2.0 * phi1 - phi[0];
          phi[1] = -phi[1];
        }
        if (phi_old <= phi2 && phi[0] > phi2) {
          phi[0] = 2.0 * phi2 - phi[0];
          phi[1] = -phi[1];
        }
      #elif MODEL==3
        NFLOOP{
          phi[2*nf] += phiamp * sqrt_dN * dwlist[nf-1][i];
        }
      #else
        phi[0] += phiamp * dw * sqrt_dN;
        phi[0] += phiamp * bias * Bias * GaussianFactor * dN;
      #endif

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
    }

    animationcount++;
    if(sanimation && animationcount>aninum-1){
      animation(phievol, seed, animationstep);
      animationstep++;
      animationcount=0;
    }

    N += dN;
    std::cout << "\rLatticeSimulation   : " << std::setw(1) << int(100.*(n+1)/(double)(endtime)) << "%" << std::flush;
  }
  std::cout << std::endl;
}


void evolutionNoise(int seed, int averagetime) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NLnoiseAll; i++) {
    double N = totalstep*dN;
    std::seed_seq seq{seed,averagetime,i};
    std::mt19937 engine_av(seq);
    std::normal_distribution<> dist_av(0., 1.);

    state_type phi = phievol[i];

    // stepper
    boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;
    
    while (EoN(phi)>0) {
      #if MODEL==3
        double psiamp = sqrt(calPpsi(phi));
        stepper_noise.do_step(dphidN, phi, N, dN);
        N += dN;
        NFLOOP{
          double dw = dist_av(engine_av);
          phi[2*nf] += psiamp * dw * sqrt_dN;
        }
      #elif MODEL==2
        double calPphival = calPphi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double calPpival = calPpi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double crosscor = RecalPphipi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double phiamp = sqrt(calPphival);
        
        double Nstep = N;
        for (int dn=0;dn<(int)divdN;dn++) {
          stepper_noise.do_step(dphidN, phi, Nstep, dN/divdN);
          Nstep += dN/divdN;
        }
        N += dN;
        
        double dw = dist_av(engine_av);
        double dwpi = dist_av(engine_av);

        double b2 = crosscor*crosscor;
        double combi1 = sqrt(4.*b2 + pw2(calPphival-calPpival));

        double sqrtlam1 = sqrt(0.5*(calPphival+calPpival + combi1));
        double eig2 = 0.5*(calPphival+calPpival - combi1);
        if(eig2 < 1e-14*(calPphival+calPpival)) eig2 = 0.;

        double sqrtlam2 = sqrt(eig2);

        double denomplus = sqrt(4.*b2 + pw2(calPphival-calPpival + combi1));
        double denomminus = sqrt(4.*b2 + pw2(calPphival-calPpival - combi1));

        double vplus1 = (calPphival-calPpival + combi1)/denomplus;
        double vplus2 = 2.*crosscor/denomplus;
        double vminus1 = (calPphival-calPpival - combi1)/denomminus;
        double vminus2 = 2.*crosscor/denomminus;

        phi[0] += (sqrtlam1*vplus1*dw + sqrtlam2*vminus1*dwpi) * sqrt_dN;
        phi[1] += (sqrtlam1*vplus2*dw + sqrtlam2*vminus2*dwpi) * sqrt_dN;
      #elif MODEL==1
        double phiamp = sqrt(calPphi(N,phi,N0list[i],brokenlist[i]));
        double piamp = sqrt(calPpi(N,phi,N0list[i],brokenlist[i]));
        double crosscor = RecalPphipi(N,phi,N0list[i],brokenlist[i]);
      #else
        double phiamp = sqrt(calPphi(phi));
        
        stepper_noise.do_step(dphidN, phi, N, dN);
        N += dN;
        
        double dw = dist(engine_av);
        phi[0] += phiamp * dw * sqrt_dN;
      #endif

      if(strajectory && i==0 && superH) {
        save_trajectory(phi, N+dN*totalstep);
      }
    }

    phievol[i] = phi;
    Nnoise[i] = N - dN*totalstep;
  }
}


void dNmap(int InterpolatingNo) {

  // guided (not static) is intentional here: each point runs its own
  // dense-output zero-crossing search below, and the number of steps to
  // converge genuinely differs per point -- unlike the drift/noise loops in
  // evolution(), this one has real per-iteration load imbalance to balance.
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
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

    int countstep = 0;
    int caouttrj=0;
    while (true){
      stepper.do_step(dphidN);
      phi = stepper.current_state();
      N = stepper.current_time();

      caouttrj++;
      if(strajectory && i==0 && caouttrj>1000){
        double Nd = N;
        save_trajectory(phi, Nd);
        caouttrj=0;
      }
      
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
      countstep++;
      if (countstep>1e4) {
        std::cout << i << " N = " << N << std::endl;
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
