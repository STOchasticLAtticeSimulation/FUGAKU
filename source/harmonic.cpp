// g++ -std=c++14 -I/opt/homebrew/opt/boost/include harmonic.cpp -o harmonic

#include <boost/numeric/odeint.hpp>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>

std::random_device seed;
std::mt19937 engine(seed());
std::normal_distribution<> dist(0., 1.);

// type
typedef std::vector<double> state_type;

// EoM
void harmonic_oscillator(const state_type &x, state_type &dxdt, const double t) {
  const double omega = 1.;
  dxdt[0] = x[1];
  dxdt[1] = -omega * omega * x[0];
}

// Observers
void output_to_file(state_type &x, double t, std::ofstream &outfile) {
  outfile << t << " " << x[0] << " " << x[1] << std::endl;
}

void AddNoise(state_type &x, double t) {
  x[0] += 0.0003*dist(engine);
}


int main() {
  std::ofstream outfile("harmonic_oscillator.dat");

  state_type x(2);
  x[0] = 1.0;
  x[1] = 0.0;

  // time
	double ti = 0.0;
  double tnoise = 2.0;
  double tend = 10.0;
  double t = ti;
  double dt = 0.01;

  // stepper
	boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;
	// auto stepperN = make_dense_output(1.0e-7, 1.0e-7, stepper_nosee());

	auto observer = [&outfile](state_type &x, double t) {
  	output_to_file(x, t, outfile);
		AddNoise(x,t);
  };
	
	integrate_const(stepper_noise, harmonic_oscillator, x, ti, tnoise, dt, observer);

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
		stepper.do_step(harmonic_oscillator);
		x = stepper.current_state();
		t = stepper.current_time();
		if (xrev[0]*x[0]<0){
			double prec_x = 10.;
			double prec = 1.e-7;
			double tl = t - stepper.current_time_step();
			double tr = t;
			double tmid = tl + (tr - tl)/2.;
			while (prec_x>prec){
				stepper.calc_state(tl, xbefore);
				stepper.calc_state(tr, xafter);
				tmid = (tl*xafter[0]-tr*xbefore[0]) / (xafter[0]-xbefore[0]);
				stepper.calc_state(tmid, xmid);
				if (xmid[0]>0) tr = tmid;
				else tl = tmid;
				prec_x = std::fabs(xmid[0]);
			}

			stepper.calc_state(tmid, xmid);
			outfile << tmid << " " << xmid[0] << " " << xmid[1] << std::endl;
			break;
		}

		outfile << t << " " << x[0] << " " << x[1] << std::endl;	
		// std::cout << tafter << " " << xafter[0] << std::endl;
		}

    return 0;
}



// #include <iostream>
// #include <random>
// #include <vector>
// #include <fstream>
// #include <boost/numeric/odeint.hpp>


// std::random_device seed;
// std::mt19937 engine(seed());
// std::normal_distribution<> dist(0., 1.);

// using namespace std;
// using namespace boost::numeric::odeint;

// typedef vector<double> state_type;

// void harmonic_oscillator(const state_type &x, state_type &dxdt, double t) {
//     const double omega = 1.0;
//     dxdt[0] = x[1];
//     dxdt[1] = -omega * omega * x[0];
// }

// void output_to_file(state_type &x, double t, ofstream &outfile) {
//   outfile << t << " " << x[0] << " " << x[1] << endl;
// }

// void AddNoise(state_type &x, double t) {
//   if(t<5) x[0] += 0.1*dist(engine);
// }


// int main() {
//     state_type x(2);
//     x[0] = 1.0;
//     x[1] = 0.0;

//     double ti = 0.0;
//     double t_end = 10.0;
//     double dt = 0.2;

//     // stepper
// 		// runge_kutta4<state_type> stepper;
// 		typedef runge_kutta_dopri5<state_type>base_stepper_type;
// 		auto stepper = make_dense_output( 1.0e-10 , 1.0e-6 , base_stepper_type());


//     // 出力ファイルのオープン
//     ofstream outfile("harmonic_oscillator.dat");

// 		auto observer = [&outfile](state_type &x, double t) {
//       output_to_file(x, t, outfile);
// 			// AddNoise(x,t);
//     };

//     integrate_const(stepper, harmonic_oscillator, x, ti, t_end, dt, observer);

//     outfile.close();

//     return 0;
// }
