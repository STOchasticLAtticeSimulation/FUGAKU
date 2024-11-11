#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <complex>
#include <algorithm>
#include <sys/time.h>

#include "parameters.hpp"
#include "fft.hpp"
#include "vec_op.hpp"


const std::string MODELdir = "../Starobinsky/";

const std::string Nfileprefix = sdatadir + "/Nmap_";
const std::string powfileprefix = sdatadir + "/power_";
const std::string cmpfileprefix = sdatadir + "/compaction_";
const std::string prbfileprefix = sdatadir + "/probabilities";
const std::string powsfileprefix = sdatadir + "/powers";
std::ifstream Nfile;
std::ofstream powfile, cmpfile, prbfile, powsfile;

// useful macro
#define LOOP for(int i = 0; i < NLnoise; i++) for(int j = 0; j < NLnoise; j++) for(int k = 0; k < NLnoise; k++)


bool compareByFirstColumn(const std::vector<double>& a, const std::vector<double>& b);
void spectrum_div(std::vector<std::vector<double>> Ndata, int noisefiledirNo);
void compaction_div(std::vector<std::vector<double>> Ndata, int noisefiledirNo);


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
  std::cout << MODELdir + Nfileprefix + std::to_string(NLnoise) + "_" + std::to_string(noisefiledirNo) << std::endl;
  Nfile.open(MODELdir + Nfileprefix + std::to_string(NLnoise) + "_" + std::to_string(noisefiledirNo) + std::string(".dat"));
  std::vector<std::vector<double>> Ndata;
  if (Nfile.fail()) std::cout << "Nfile fail" << "\n";

  // データの読み込み
  double first, second;
  while (Nfile >> first >> second) {
    Ndata.push_back({first, second});
  }
  Nfile.close();

  // データのソート
  std::sort(Ndata.begin(), Ndata.end(), compareByFirstColumn);

  if(spower) spectrum_div(Ndata,noisefiledirNo);
  if(scompaction) compaction_div(Ndata,noisefiledirNo);


  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  std::cout << std::endl;
  // -------------------------------------
}


// 比較関数: 1列目（整数部分）を基準にソート
bool compareByFirstColumn(const std::vector<double>& a, const std::vector<double>& b) {
  return a[0] < b[0];
}


// Calculate power spectrum
void spectrum_div(std::vector<std::vector<double>> Ndata, int noisefiledirNo) {
  // Nfile.open(MODELdir + Nfileprefix + std::to_string(NLnoise) + "_" + std::to_string(noisefiledirNo) + std::string(".dat"));
  // std::vector<std::vector<double>> Ndata;
  // if (Nfile.fail()) std::cout << "Nfile fail" << "\n";

  // // データの読み込み
  // double first, second;
  // while (Nfile >> first >> second) {
  //   Ndata.push_back({first, second});
  // }
  // Nfile.close();

  // // データのソート
  // std::sort(Ndata.begin(), Ndata.end(), compareByFirstColumn);

  powfile.open(MODELdir + powfileprefix + std::to_string(NLnoise) + std::string("_") + std::to_string(noisefiledirNo) + std::string(".dat"));
  powfile << std::setprecision(10);
  std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D = std::vector<std::vector<std::vector<std::complex<double>>>>(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise, 0)));

  for (int i=0; i<NLnoise*NLnoise*NLnoise; i++) {
    int x=i/NLnoise/NLnoise ,y=(i%(NLnoise*NLnoise))/NLnoise, z=i%NLnoise;
    Nmap3D[x][y][z] = Ndata[i][1];
  }
  std::vector<std::vector<std::vector<std::complex<double>>>> Nk=fft(Nmap3D);
  
  powsfile.open(MODELdir + powsfileprefix + std::string(".dat"), std::ios::app);
  powsfile << std::setprecision(10);
  int imax = ceil(log(NLnoise/2)/dlogn);
  std::vector<double> disc_power(imax, 0);

  LOOP{
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
void compaction_div(std::vector<std::vector<double>> Ndata, int noisefiledirNo) {
  // Nfile.open(MODELdir + Nfileprefix + std::to_string(NLnoise) + "_" + std::to_string(noisefiledirNo) + std::string(".dat"));
  // std::vector<std::vector<double>> Ndata;
  // std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D = std::vector<std::vector<std::vector<std::complex<double>>>>(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise, 0)));
  // if (Nfile.fail()) std::cout << "Nfile fail" << "\n";
  // std::cout << MODELdir + Nfileprefix + std::to_string(NLnoise) + "_" + std::to_string(noisefiledirNo) << std::endl;

  // // データの読み込み
  // double first, second;
  // while (Nfile >> first >> second) {
  //   Ndata.push_back({first, second});
  // }
  // Nfile.close();

  // // データのソート
  // std::sort(Ndata.begin(), Ndata.end(), compareByFirstColumn);

  prbfile.open(MODELdir + prbfileprefix + std::string(".dat"), std::ios::app);
  cmpfile.open(MODELdir + cmpfileprefix + std::to_string(NLnoise) + std::string("_") + std::to_string(noisefiledirNo) + std::string(".dat"));
  prbfile << std::setprecision(10);
  cmpfile << std::setprecision(10);
  
  double Naverage = 0;
  int dr = 1;
  for (size_t n = 0; n < Ndata.size(); n++) {
    Naverage += Ndata[n][1];
  }
  Naverage /= NLnoise*NLnoise*NLnoise;

  // zeta map
  for (size_t n = 0; n < Ndata.size(); n++) {
    Ndata[n][1] -= Naverage;
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
        zetar[1][ri]+=Ndata[i][1];
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
