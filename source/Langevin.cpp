// g++ -std=c++14 -O2 -I/opt/homebrew/opt/boost/include Langevin.cpp -o Langevin

#include <boost/numeric/odeint.hpp>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>

const double mm = 0.0211;
int output_step = 0;

// Noise term
std::random_device seed;
std::mt19937 engine(seed());
std::normal_distribution<> dist(0., 1.);

// type
typedef std::vector<double> state_type;


// Potential
double VV(double phi) {
  return mm*mm*phi*phi/2.;
}

// Derivative
double Vp(double phi) {
  return mm*mm*phi;
}

double hubble(double phi, double pi) {
  return sqrt((pi*pi/2. + VV(phi))/3.);
}

double epm1(double phi, double pi) {
  double HH = hubble(phi,pi);
  return pi*pi/2./HH/HH - 1.;
}


// EoM
void dphidN(const state_type &x, state_type &dxdt, const double t) {
  double xx = x[0]; // phi
  double pp = x[1]; // pi
  double HH = hubble(xx,pp);

  dxdt[0] = pp/HH;
  dxdt[1] = -3*pp - Vp(xx)/HH;
}



// Observers
void output_to_file(state_type &x, double t, std::ofstream &outfile) {
	if(output_step%10==0){
		outfile << t << " " << x[0] << " " << x[1] << " " << hubble(x[0],x[1]) << std::endl;
	}
	output_step++;
}

void AddNoise(state_type &x, double t) {
  x[0] += 0.0*dist(engine);
}



int main() {
  std::ofstream outfile("Langevin.dat");

  state_type x(2);
  x[0] = 11.0;
  x[1] = 0.;//-sqrt(2./3)*mm;

  // time
  double ti = 0.0;
  double tnoise = 5.0;
  double tend = 33.0;
  double t = ti;
  double dt = 0.01;

  // stepper
	boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;

	auto observer = [&outfile](state_type &x, double t) {
  	output_to_file(x, t, outfile);
		AddNoise(x,t);
  };
	
	integrate_const(stepper_noise, dphidN, x, ti, tnoise, dt, observer);

	std::cout << "Noise end" << std::endl;

	// Find 0 crossing time
	typedef boost::numeric::odeint::runge_kutta_dopri5<state_type>base_stepper_type;
	auto stepper = make_dense_output(1.0e-7, 1.0e-7, base_stepper_type());

	stepper.initialize(x, tnoise, dt);
	state_type xrev(2);
	state_type xbefore(2);
	state_type xafter(2);
	state_type xmid(2);

	while (t<tend){
		xrev = x;
		stepper.do_step(dphidN);
		x = stepper.current_state();
		t = stepper.current_time();
		
		if (epm1(xrev[0],xrev[1])*epm1(x[0],x[1])<0){
			double prec_x = 10.;
			double prec = 1.e-7;
			double tl = t - stepper.current_time_step();
			double tr = t;
			double tmid = 0;

			while (prec_x>prec){
				stepper.calc_state(tl, xbefore);
				stepper.calc_state(tr, xafter);
				tmid = (tl*epm1(xafter[0],xafter[1])-tr*epm1(xbefore[0],xbefore[1])) / (epm1(xafter[0],xafter[1])-epm1(xbefore[0],xbefore[1]));
				stepper.calc_state(tmid, xmid);
				if (epm1(xmid[0],xmid[1])>0) tr = tmid;
				else tl = tmid;
				prec_x = std::fabs(epm1(xmid[0],xmid[1]));
			}

			stepper.calc_state(tmid, xmid);
			outfile << tmid << " " << xmid[0] << " " << xmid[1] << " " << hubble(xmid[0],xmid[1]) << std::endl;
			break;
		}

		outfile << t << " " << x[0] << " " << x[1] << " " << hubble(x[0],x[1]) << std::endl;
		}

  return 0;
}
