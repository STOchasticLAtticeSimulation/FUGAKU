#ifndef INCLUDED_model_hpp_
#define INCLUDED_model_hpp_

// type
typedef std::vector<double> state_type;

// Model parameters
const std::string model = "chaotic"; // Name of the model
const double mm = 0.0211; // Mass of the inflaton
const std::vector<double> phii{11.,-sqrt(2./3)*mm}; // Initial conditions {field,derivative}


// Potential
double VV(double phi) {
  return mm*mm*phi*phi/2.;
}

// Derivative
double Vp(double phi) {
  return mm*mm*phi;
}

double hubble(const state_type &phi) {
  return sqrt((phi[1]*phi[1]/2. + VV(phi[0]))/3.);
}

double ep(const state_type &phi) {
  double HH = hubble(phi);
  return phi[1]*phi[1]/2./HH/HH;
}

// Power spectrum of phi
double calPphi(const state_type &phi) {
  double HH = hubble(phi);
  double eps = Vp(phi[0])*Vp(phi[0])/VV(phi[0])/VV(phi[0])/2.;
  double Hi = hubble(phii);
  double NoiseNLO = (1. + eps*(6.-4.*euler_gamma-8.*log(2.))) * pow(sigma*Hi/2./HH, -4.*eps);

  return pow(HH/2./M_PI, 2) * NoiseNLO;
}

// The condition at the end of inflation
inline double EoI(const state_type &phi) {
  return 1-ep(phi);
}

// EoM
void dphidN(const state_type &phi, state_type &dxdt, const double t) {
  double xx = phi[0]; // phi
  double pp = phi[1]; // pi
  double HH = hubble(phi);

  dxdt[0] = pp/HH;
  dxdt[1] = -3*pp - Vp(xx)/HH;
}

#endif