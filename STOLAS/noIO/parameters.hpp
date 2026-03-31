#ifndef INCLUDED_parameters_hpp_
#define INCLUDED_parameters_hpp_

// Parameters of STOLAS
const double sigma = pow(2.,-4.); //0.1; // //ksigma = 2pi sigma exp(N) / L, nsigma = sigma exp(N)
constexpr double dn = 1; // Thickness of nsigma sphere shell
constexpr int NLnoise = 64; // Box size L
constexpr int totalnoiseNo = 1;//pow(2,3); // The number of chunks
constexpr double dN = 0.01*LOG2; // e-folds step
constexpr double Nprec = 1e-7; // Precision of e-foldings
constexpr double dlogn = 0.1; // Width of bin in power spectrum
constexpr double LL = 1.; // Box size L

double nsigmareset = 8.;
constexpr int aninum = 20; // dividing number for animation
constexpr bool EoI_noise = false; // Adding the noise until EoN
constexpr int MeanNumber = 1; // Number of averageing in each
constexpr int internumber = 0; // Number of zoom-in

// Outputs
constexpr bool sfield = false; // Output field
constexpr bool strajectory = false; // Output the trajectory
constexpr bool spower = true; // Output the power spectrum
constexpr bool sanimation = false; // Output the animation
constexpr bool sweight = false; // Output the weight
constexpr bool scompaction = false; // Output the compaction

// Importance sampling
constexpr double Nbias = 3.8; // Time of the bias
constexpr double dNbias = 0.1; // Variance of the bias
const double bias = 0*sqrt(dNbias); // Amplitude of the bias

// Directory name of saved data, you can change after "make clean" in your terminal.
const std::string sdatadir = "data";

#endif
