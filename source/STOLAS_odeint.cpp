#include "STOLAS.hpp"
#include "vec_op.hpp"
#include "fft.hpp"

// -- Transpose ---------------------------------
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix) {
  if (matrix.empty() || matrix[0].empty()) return {}; // 空の場合の処理

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

// -- Odeint ------------------------------------
#include <boost/numeric/odeint.hpp>

// type
typedef std::vector<double> state_type;

// Model parameters
const double H0 = 1e-5; // Hubble parameter of broken point
const double calPRIR = 8.5e-10; // Amplitude of curvature perturbation
const double Lambda = 1700; // Ratio between Ap to Am
const double Ap = sqrt(9./4/M_PI/M_PI*H0*H0*H0*H0*H0*H0/calPRIR); // Gradient of the potential at first stage
const double Am = Ap/Lambda; // Gradient of the potential at second stage
const double V0 = 3*H0*H0; // Amplitude of flat potential
const double phif = -0.0187; // The inflaton value at the end of inflation

// Potential
double VV(double phi) {
  if (phi > 0) {
    return V0 + Ap*phi;
  } else {
    return V0 + Am*phi;
  }
}

// Derivative
double Vp(double phi) {
  if (phi > 0) {
    return Ap;
  } else {
    return Am;
  }
}

double hubble(double phi, double pi) {
  return sqrt((pi*pi/2. + VV(phi))/3.);
}

// Power spectrum of phi
double calPphi(double &N, std::vector<double> &phi, double N0, bool broken) {
  if (!broken) {
    return pow(hubble(phi[0],phi[1])/2./M_PI,2);
  } else {
    double alpha = exp(N-N0);
    return pow(hubble(phi[0],phi[1])/2./M_PI,2) *
      ((1 + pow(sigma,2))*(9*pow(1 + pow(alpha,2)*pow(sigma,2),2) - 
				18*Lambda*pow(1 + pow(alpha,2)*pow(sigma,2),2) + 
				pow(Lambda,2)*(9 + 18*pow(alpha,2)*pow(sigma,2) + 9*pow(alpha,4)*pow(sigma,4) + 
					       2*pow(alpha,6)*pow(sigma,6))) + 
	    3*(-3*(1 + (-1 + 4*alpha)*pow(sigma,2) - (-4 + alpha)*pow(alpha,3)*pow(sigma,4) + 
		   pow(alpha,4)*pow(sigma,6)) + Lambda*
	       (6 + 6*(-1 + 4*alpha)*pow(sigma,2) + 2*(14 - 5*alpha)*pow(alpha,3)*pow(sigma,4) + 
		2*(5 - 2*alpha)*pow(alpha,4)*pow(sigma,6)) + 
	       pow(Lambda,2)*(-3 + (3 - 12*alpha)*pow(sigma,2) + pow(alpha,3)*(-16 + 7*alpha)*pow(sigma,4) + 
			      pow(alpha,4)*(-7 + 4*alpha)*pow(sigma,6)))*cos(2*(-1 + alpha)*sigma) + 
	    6*sigma*(-3*pow(-1 + Lambda,2) + pow(alpha,4)*(3 - 10*Lambda + 7*pow(Lambda,2))*pow(sigma,4) - 
		     3*alpha*pow(-1 + Lambda,2)*(-1 + pow(sigma,2)) - 
		     pow(alpha,3)*(3 - 7*Lambda + 4*pow(Lambda,2))*pow(sigma,2)*(-1 + pow(sigma,2)) + 
		     pow(alpha,5)*(-1 + Lambda)*Lambda*pow(sigma,4)*(-1 + pow(sigma,2)))*sin(2*sigma - 2*alpha*sigma))/
      (2.*pow(alpha,6)*pow(Lambda,2)*pow(sigma,6));
  }
}

// The power spectrum of pi
double calPpi(double &N, std::vector<double> &phi, double N0, bool broken) {
  return 0;
}

// Cross correlation of phi and pi
double RecalPphipi(double &N, std::vector<double> &phi, double N0, bool broken) {
  return 1;
}

// The condition at the end of inflation
inline double EoI(state_type &phi) {
  return phi[0] - phif;
}


// EoM
void dphidN(const state_type &x, state_type &dxdt, const double t) {
  double xx = x[0]; // phi
  double pp = x[1]; // pi
  double HH = hubble(xx,pp);

  dxdt[0] = pp/HH;
  dxdt[1] = -3*pp - Vp(xx)/HH;
}

// ----------------------------------------------


// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)


STOLAS::STOLAS(std::string Model, double DN, std::string sourcedir, int NoisefileDirNo, std::vector<double> Phii, double Bias, double NBias, double DNbias, int NoisefileNo) {

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif
  
  model = Model;
  dN = DN;
  noisefiledirNo = NoisefileDirNo;
  phii = Phii;
  bias = Bias;
  Nbias = NBias;
  dNbias = DNbias;
  noisefileNo = NoisefileNo;

  std::cout << "Noise file No. : " << noisefileNo << std::endl;
  
  noisefile.open(sourcedir + std::string("/") + noisefiledir + std::to_string(noisefiledirNo) + noisefilenamediv + std::to_string(noisefileNo) + std::string(".bin"), std::ios::binary);
  noisefilefail = noisefile.fail();
  biasfile.open(sourcedir + std::string("/") + biasfilenamediv + std::to_string(noisefileNo) + std::string(".bin"), std::ios::binary);
  biasfilefail = biasfile.fail();
  
  if (!noisefile.fail() && !biasfile.fail()) {
    std::cout << "model : " << model << std::endl;

    int totalstep = ceil(log((NLnoise/2-1)/sigma)/dN);
    while (true) {
      std::vector<double> vv(NL);
      noisefile.read(reinterpret_cast<char*>(vv.data()), NL * sizeof(double));

      // EOFに到達して完全に読み込まれなかった場合を確認
      if (noisefile.gcount() != NL * sizeof(double)) {
        break; // 全行の読み込み完了または不完全な行をスキップ
      }
      noisedataTr.push_back(vv);
    }
    noisedata = transpose(noisedataTr);
    noisedataTr.clear();

    while (true) {
      std::vector<double> vv(NL);
      biasfile.read(reinterpret_cast<char*>(vv.data()), NL * sizeof(double));

      // EOFに到達して完全に読み込まれなかった場合を確認
      if (biasfile.gcount() != NL * sizeof(double)) {
        break; // 全行の読み込み完了または不完全な行をスキップ
      }
      biasdataTr.push_back(vv);
    }
    biasdata = transpose(biasdataTr);
    biasdataTr.clear();


    std::cout << "Noise/Bias data imported. Box size is " << NL << "." << std::endl;
    Nfile.open(Noiseprefix + std::to_string(NLnoise) + std::string("_") + std::to_string(NoisefileDirNo) + std::string(".dat"), std::ios::app);
    logwfile.open(logwfileprefix + std::to_string(NLnoise) + std::string(".dat"), std::ios::app);

    Hdata = std::vector<std::vector<double>>(noisedata[0].size(), std::vector<double>(NL,0));
    pidata = std::vector<std::vector<double>>(noisedata[0].size(), std::vector<double>(NL,0));
    Ndata = std::vector<double>(NL,0);
  }
}


bool STOLAS::checknoisefile() {
  return !noisefilefail;
}

bool STOLAS::checkbiasfile() {
  return !biasfilefail;
}

bool STOLAS::noisebiassize() {
  return (noisedata.size() == biasdata.size());
}

bool STOLAS::Nfilefail() {
  return Nfile.fail();
}


void STOLAS::dNmap(int NoiseNo) {
  Nfile << std::setprecision(10);
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
        Hdata[n][i] = pow(hubble(phi[0],phi[1]),2);
        pidata[n][i] = phi[1]*phi[1];
      }
      double phiamp = sqrt(calPphi(N,phi,N0,broken));
      double piamp = sqrt(calPpi(N,phi,N0,broken));
      double crosscor = RecalPphipi(N,phi,N0,broken);
      
      double dw = noisedata[i][n];
      double Bias = biasdata[i][n];

      stepper_noise.do_step(dphidN, phi, N, dN);
      N += dN;
      
      phi[0] += phiamp * dw * sqrt(dN);

      if (crosscor > 0) {
        phi[1] += piamp * dw * sqrt(dN);
      } else {
        phi[1] -= piamp * dw * sqrt(dN);
      }

      double GaussianFactor = 1./dNbias/sqrt(2*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
      phi[0] += phiamp * bias * Bias * GaussianFactor * dN;

      if (crosscor > 0) {
        phi[1] += piamp * bias * Bias * GaussianFactor * dN;
      } else {
        phi[1] -= piamp * bias * Bias * GaussianFactor * dN;
      }

      if (!broken && phi[0] < 0) {
        N0 = N;
        broken = true;
      }
      
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
        double prec = 1.e-7;
        double Nl = N;
        double Nr = N - stepper.current_time_step();
        double Nmid = 0;

        while (precphi>prec){
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
      Ndata[i] = N;
      Nfile << i + NL*NoiseNo << ' ' << N << std::endl;
      complete++;
      std::cout << "\rLatticeSimulation : " << std::setw(3) << 100*complete/NL << "%" << std::flush;
    }
  }
  std::cout << std::endl;

}



// calculation of weight
void STOLAS::weight() {
  double logw = 0.;
  for (size_t n=0; n<noisedata[0].size(); n++) {
    double N = n*dN;
    double Bias = bias/dNbias/sqrt(2*M_PI)*exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
    logw -= Bias*noisedata[0][n]*sqrt(dN) + (Bias*Bias*dN)/2;
  }
  logwfile << noisefileNo << ' ' << logw << std::endl;
}

// Export animation
void STOLAS::animation() {
  Hfile.open(Hfileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
  pifile.open(pifileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
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
}
