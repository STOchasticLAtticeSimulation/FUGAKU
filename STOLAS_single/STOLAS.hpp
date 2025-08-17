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
#include "src/fft.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)
#define LOOPLONG for(int i = 0; i < NLnoise; i++) for(int j = 0; j < NLnoise; j++) for(int k = 0; k < NLnoise; k++)


const std::string Nfileprefix = sdatadir + "/" + model + "/Nmap_";
const std::string Hfileprefix = sdatadir + "/" + model + "/animation/H_";
const std::string pifileprefix = sdatadir + "/" + model + "/animation/pi_";
const std::string powfileprefix = sdatadir + "/" + model + "/power_";
const std::string cmpfileprefix = sdatadir + "/" + model + "/compaction_";
const std::string prbfileprefix = sdatadir + "/" + model + "/probabilities";
const std::string powsfileprefix = sdatadir + "/" + model + "/powers";
const std::string logwfileprefix = sdatadir + "/" + model + "/logw_";
bool noisefilefail, biasfilefail, Nfilefail;


int Nn, noisefiledirNo, noisefileNo;
std::ofstream Nfile, Hfile, pifile, powfile, cmpfile, prbfile, powsfile, logwfile;
std::vector<std::vector<double>> noisedata, noisedataTr, biasdata, biasdataTr, Hdata, pidata;
std::vector<double> Ndata;


bool checknoisefile();
bool checkbiasfile();
bool noisebiassize();
bool checkNfilefail();
 
void dNmap(int noisefileNo);
void animation(int NoiseNo);
void weight();
void spectrum(std::vector<double> Ndata, int noisefiledirNo);
void compaction(std::vector<double> Ndata, int noisefiledirNo);

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix);
bool compareByFirstColumn(const std::vector<double>& a, const std::vector<double>& b);


void STOLAS(int NoisefileDirNo, int NoisefileNo);



void STOLAS(int NoisefileDirNo, int NoisefileNo) {
  std::ifstream noisefile, biasfile;

  std::cout << "Noise file No. : " << NoisefileNo << std::endl;
  
  noisefile.open(noisefiledir + std::to_string(NoisefileDirNo) + noisefilenamediv + std::to_string(NoisefileNo) + std::string(".bin"), std::ios::binary);
  noisefilefail = noisefile.fail();
  biasfile.open(biasfilenamediv + std::to_string(NoisefileNo) + std::string(".bin"), std::ios::binary);
  biasfilefail = biasfile.fail();
  
  if (!noisefile.fail() && !biasfile.fail()) {
    std::cout << "model : " << model << std::endl;

    while (true) {
      std::vector<double> vv(NL);
      noisefile.read(reinterpret_cast<char*>(vv.data()), NL * sizeof(double));

      if (noisefile.gcount() != NL * sizeof(double)) {
        break;
      }
      noisedataTr.push_back(vv);
    }
    noisedata = transpose(noisedataTr);
    noisedataTr.clear();

    while (true) {
      std::vector<double> vv(NL);
      biasfile.read(reinterpret_cast<char*>(vv.data()), NL * sizeof(double));

      if (biasfile.gcount() != NL * sizeof(double)) {
        break;
      }
      biasdataTr.push_back(vv);
    }
    biasdata = transpose(biasdataTr);
    biasdataTr.clear();


    std::cout << "Noise/Bias data imported. Box size is " << NL << "." << std::endl;
  }
}


void dNmap(int NoiseNo) {
  int complete = 0;
  int totalstep = ceil(log((NLnoise/2-1)/sigma)/dN);
  double Nnoise = (double)totalstep*dN;
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL; i++) {
    double N=0;
    int numstep = 0;
    #if MODEL==1
      double N0;
      bool broken = false;
    #endif
    state_type phi = phii;

    // stepper
    boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;

    for (size_t n=0; n<totalstep; n++) {
      if (sanimation) {
        Hdata[n][i + NL*NoiseNo] = pow(hubble(phi),2.);
        pidata[n][i + NL*NoiseNo] = phi[1]*phi[1];
      }

      #if MODEL==0
        double phiamp = sqrt(calPphi(phi));
      #elif MODEL==1
        double phiamp = sqrt(calPphi(N,phi,N0,broken));
        double piamp = sqrt(calPpi(N,phi,N0,broken));
        double crosscor = RecalPphipi(N,phi,N0,broken);
      #endif

      double dw = noisedata[i][n];
      double Bias = biasdata[i][n];
      double GaussianFactor = 1./dNbias/sqrt(2*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);

      stepper_noise.do_step(dphidN, phi, N, dN);
      N += dN;
      
      phi[0] += phiamp * dw * sqrt(dN);
      phi[0] += phiamp * bias * Bias * GaussianFactor * dN;

      #if MODEL==1
        if (crosscor > 0) {
          phi[1] += piamp * dw * sqrt(dN);
          phi[1] += piamp * bias * Bias * GaussianFactor * dN;
        } else {
          phi[1] -= piamp * dw * sqrt(dN);
          phi[1] -= piamp * bias * Bias * GaussianFactor * dN;
        }

        if (!broken && phi[0] < 0) {
          N0 = N;
          broken = true;
        }
      #endif
      
    }

    // Find zero crossing time
    typedef boost::numeric::odeint::runge_kutta_dopri5<state_type>base_stepper_type;
    auto stepper = make_dense_output(1.0e-7, 1.0e-7, base_stepper_type());

    stepper.initialize(phi, N, dN);
    state_type phil(2);
    state_type phir(2);
    state_type phim(2);

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
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      Nfile << i + NL*NoiseNo << ' ' << N << std::endl;
      if (spower || scompaction) Ndata[i + NL*NoiseNo] = N;
      complete++;
      std::cout << "\rLatticeSimulation : " << std::setw(3) << 100*complete/NL << "%" << std::flush;
    }
  }
  std::cout << std::endl;

}


bool checknoisefile() {
  return !noisefilefail;
}

bool checkbiasfile() {
  return !biasfilefail;
}

bool noisebiassize() {
  return (noisedata.size() == biasdata.size());
}

bool checkNfilefail() {
  return Nfilefail;
}


// calculation of weight
void weight() {
  double logw = 0.;
  for (size_t n=0; n<noisedata[0].size(); n++) {
    double N = n*dN;
    double Bias = bias/dNbias/sqrt(2*M_PI)*exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
    logw -= Bias*noisedata[0][n]*sqrt(dN) + (Bias*Bias*dN)/2;
  }
  logwfile << noisefileNo << ' ' << logw << std::endl;
}


// Export animation
void animation() {
  Hfile.open(Hfileprefix + std::to_string(NLnoise) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
  pifile.open(pifileprefix + std::to_string(NLnoise) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
  Hfile << std::setprecision(14);
  pifile << std::setprecision(14);

  for (size_t n=0; n<Hdata.size(); n++) {
    for (size_t i=0; i<Hdata[n].size(); i++) {
      Hfile << Hdata[n][i] << ' ';
      pifile << pidata[n][i] << ' ';
    }
    Hfile << std::endl;
    pifile << std::endl;
    std::cout << "\rAnimeDataExporting : " << std::setw(3) << 100*(n+1)/Hdata.size() << "%" << std::flush;
  }
  std::cout << std::endl;
  Hfile.close();
  pifile.close();
}


std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix) {
  if (matrix.empty() || matrix[0].empty()) return {};

  size_t rows = matrix.size();
  size_t cols = matrix[0].size();
  std::vector<std::vector<double>> result(cols, std::vector<double>(rows));

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      result[j][i] = matrix[i][j];
    }
  }
  return result;
}


bool compareByFirstColumn(const std::vector<double>& a, const std::vector<double>& b) {
  return a[0] < b[0];
}


// Calculate power spectrum
void spectrum(std::vector<double> Ndata, int noisefiledirNo) {
  powfile.open(powfileprefix + std::to_string(NLnoise) + std::string("_") + std::to_string(noisefiledirNo) + std::string(".dat"));
  powfile << std::setprecision(10);
  std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D = std::vector<std::vector<std::vector<std::complex<double>>>>(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise, 0)));

  for (int i=0; i<NLnoise*NLnoise*NLnoise; i++) {
    int x=i/NLnoise/NLnoise ,y=(i%(NLnoise*NLnoise))/NLnoise, z=i%NLnoise;
    Nmap3D[x][y][z] = Ndata[i];
  }
  std::vector<std::vector<std::vector<std::complex<double>>>> Nk=fft_fftw(Nmap3D);
  
  powsfile.open(powsfileprefix + std::string(".dat"), std::ios::app);
  powsfile << std::setprecision(10);
  int imax = ceil(log(NLnoise/2)/dlogn);
  std::vector<double> disc_power(imax, 0);

  LOOPLONG{
    int nxt, nyt, nzt; // shifted index
    if (i<=NLnoise/2) {
      nxt = i;
    } else {
      nxt = i-NLnoise;
    }

    if (j<=NLnoise/2) {
      nyt = j;
    } else {
      nyt = j-NLnoise;
    }

    if (k<=NLnoise/2) {
      nzt = k;
    } else {
      nzt = k-NLnoise;
    }
    
    double rk=nxt*nxt+nyt*nyt+nzt*nzt;
    powfile << sqrt(rk) << " " << norm(Nk[i][j][k])/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise << std::endl;

    double LogNk = log(sqrt(rk));
    double calPk = norm(Nk[i][j][k])/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise;
    for (size_t ii = 0; ii < imax; ii++) {
      if (std::abs(dlogn*ii-LogNk)<=dlogn/2.) {
        disc_power[ii] += calPk/dlogn;
        break;
      }
    }
  }
  powsfile << noisefiledirNo << " ";
  for (size_t ii = 0; ii < imax; ii++) {
    powsfile << disc_power[ii] << " " ;
  }
  powsfile << std::endl;
  std::cout << "ExportPowerSpectrum" << std::endl;
}


// Calculate compaction function
void compaction(std::vector<double> Ndata, int noisefiledirNo) {
  prbfile.open(prbfileprefix + std::string(".dat"), std::ios::app);
  cmpfile.open(cmpfileprefix + std::to_string(NLnoise) + std::string("_") + std::to_string(noisefiledirNo) + std::string(".dat"));
  prbfile << std::setprecision(10);
  cmpfile << std::setprecision(10);
  
  double Naverage = 0;
  int dr = 1;
  for (size_t n = 0; n < Ndata.size(); n++) {
    Naverage += Ndata[n];
  }
  Naverage /= NLnoise*NLnoise*NLnoise;

  // zeta map
  for (size_t n = 0; n < Ndata.size(); n++) {
    Ndata[n] -= Naverage;
  }

  // radial profile
  std::vector<std::vector<double>> zetar(2, std::vector<double>(NLnoise/2,0));
  for (size_t i=0; i<NLnoise*NLnoise*NLnoise; i++) {
    int nx=i/NLnoise/NLnoise ,ny=(i%(NLnoise*NLnoise))/NLnoise, nz=i%NLnoise;

    // centering
    if (nx<=NLnoise/2) {
      nx = nx;
    }
    else {
      nx = nx-NLnoise;
    }
    if (ny<=NLnoise/2) {
      ny = ny;
    }
    else {
      ny = ny-NLnoise;
    }
    if (nz<=NLnoise/2) {
      nz = nz;
    }
    else {
      nz = nz-NLnoise;
    }

    for (size_t ri=0; ri<NLnoise/2; ri++) {
      double norm = std::abs(sqrt(nx*nx+ny*ny+nz*nz)-ri);
      if (norm<=dr/2.) {
        zetar[0][ri]++;
        zetar[1][ri]+=Ndata[i];
        break;
      }
    }
  }
  for (size_t ri=0; ri<NLnoise/2; ri++) {
    zetar[1][ri] /= zetar[0][ri]; // average
  }

  // derivative zeta
  std::vector<double> dzetar(NLnoise/2,0);
  for (size_t ri=1; ri<NLnoise/2-1; ri++) {
    dzetar[ri] = (zetar[1][ri+1] - zetar[1][ri-1])/(2.*dr);
  }

  // compaction function
  double CompactionMax=0, CompactionInt=0, rmax=0, Rmax=0, IntTemp=0;
  bool Cnegative = false;
  for (size_t ri=0; ri<NLnoise/2; ri++) {
    double CompactionTemp = 2./3.*(1. - pow(1 + ri*dzetar[ri], 2));
    IntTemp += ri*ri*CompactionTemp*exp(3.*zetar[1][ri])*(1 + ri*dzetar[ri]);

    if (!Cnegative) {
      if (CompactionTemp < -0.2) {
	Cnegative = true;
      }
      
      if (CompactionMax<CompactionTemp) {
	CompactionMax = CompactionTemp;
	rmax = ri;
	Rmax = exp(zetar[1][ri])*ri;
	CompactionInt += IntTemp;
	IntTemp = 0;
      }
    }
    
    cmpfile << ri << ' ' << CompactionTemp << std::endl;
  }
  CompactionInt /= pow(Rmax, 3)/3.;

  prbfile << noisefiledirNo //<< ' ' << logw 
  << ' ' << CompactionInt << ' ' << CompactionMax << ' ' << Rmax << ' ' << rmax << std::endl;
  std::cout << "ExportCompactionFunction" << std::endl;
}

#endif
