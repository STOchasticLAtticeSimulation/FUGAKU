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

#include "parameters.hpp"

#define euler_gamma 0.57721566490153286061

#ifdef _OPENMP
#include <omp.h>
#endif

class STOLAS
{
protected:
  // const std::string Nfileprefix = sdatadir + "/Nmap_";
  const std::string Noiseprefix = sdatadir + "/Nmap_";
  const std::string Nfileprefix = "/Nmap_";
  const std::string Hfileprefix = sdatadir + "/H_";
  const std::string pifileprefix = sdatadir + "/pi_";
  const std::string powfileprefix = sdatadir + "/power_";
  const std::string cmpfileprefix = sdatadir + "/compaction_";
  const std::string prbfileprefix = sdatadir + "/probabilities";
  const std::string powsfileprefix = sdatadir + "/powers";
  const std::string logwfileprefix = sdatadir + "/logw_";
  bool noisefilefail, biasfilefail;

  std::string model;
  int Nn, noisefiledirNo, noisefileNo;
  double dN, bias, Nbias, dNbias;
  std::ifstream noisefile, biasfile;
  std::ofstream Nfile, Hfile, pifile, powfile, cmpfile, prbfile, powsfile, logwfile;
  std::vector<double> phii;
  std::vector<std::vector<double>> noisedata, biasdata, Hdata, pidata;
  std::vector<double> Ndata;
  std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D;

public:
  STOLAS(){}
  STOLAS(std::string Model, double DN, std::string sourcedir, int NoisefileDirNo, std::vector<double> Phii, double Bias, double NBias, double DNbias, int NoisefileNo);

  bool checknoisefile();
  bool checkbiasfile();
  bool noisebiassize();
  bool Nfilefail();
  
  double VV(double phi);
  double Vp(double phi);
  bool EoI(std::vector<double> &phi);
 
  void dNmap(int noisefileNo);
  void animation();
  void spectrum();
  void compaction();
  void weight(int noisefileNo);

  double ep(double phi, double pi);
  double hubble(double phi, double pi);

  std::vector<double> dphidN(double N, std::vector<double> phi);

  void RK4(double &t, std::vector<double> &x, double dt);

  #if MODEL==0
    double calPphi(std::vector<double> &phi);
    void RK4Mbias(double &N, std::vector<double> &phi, double dN, double dw, double Bias);
  #elif MODEL==1
    double calPphi(double &N, std::vector<double> &phi, double N0, bool broken);
    double calPpi(double &N, std::vector<double> &phi, double N0, bool broken);
    double RecalPphipi(double &N, std::vector<double> &phi, double N0, bool broken);
    void RK4Mbias(double &N, std::vector<double> &phi, double dN, double dw, double Bias, double N0, bool broken);
  #endif
};

#endif
