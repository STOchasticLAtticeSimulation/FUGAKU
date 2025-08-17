#ifndef INCLUDED_model_hpp_
#define INCLUDED_model_hpp_

// type
typedef std::vector<double> state_type;

// Model parameters
const std::string model = "Starobinsky"; // Name of the model
const double H0 = 1e-5; // Hubble parameter of broken point
const double calPRIR = 8.5e-10; // Amplitude of curvature perturbation
const double Lambda = 1700; // Ratio between Ap to Am
const double Ap = sqrt(9./4/M_PI/M_PI*H0*H0*H0*H0*H0*H0/calPRIR); // Gradient of the potential at first stage
const double Am = Ap/Lambda; // Gradient of the potential at second stage
const double V0 = 3*H0*H0; // Amplitude of flat potential
const double phif = -0.0187; // The inflaton value at the end of inflation
const std::vector<double> phii{0.0193,-5.45e-7}; // Initial conditions {field,derivative}

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

double hubble(const state_type &phi) {
  return sqrt((phi[1]*phi[1]/2. + VV(phi[0]))/3.);
}

// Power spectrum of phi
double calPphi(double &N, const state_type &phi, double N0, bool broken) {
  if (!broken) {
    return pow(hubble(phi)/2./M_PI,2);
  } else {
    double alpha = exp(N-N0);
    return pow(hubble(phi)/2./M_PI,2) *
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
double calPpi(double &N, const state_type &phi, double N0, bool broken) {
  return 0;
}

// Cross correlation of phi and pi
double RecalPphipi(double &N, const state_type &phi, double N0, bool broken) {
  return 1;
}

// The condition at the end of inflation
inline double EoI(const state_type &phi) {
  return phi[0] - phif;
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