#ifndef INCLUDED_parameters_hpp_
#define INCLUDED_parameters_hpp_

// Parameters of STOLAS
const double sigma = pow(2.,-3.); // ksigma = 2pi sigma exp(N) / L, nsigma = sigma exp(N)
const double dn = 1; // Thickness of nsigma sphere shell
constexpr int NLnoise = 64; // Box size L
const int totalnoiseNo = pow(2,3); // The number of chunks
const double dN = 0.01/log2(exp(1)); // e-folds step
const double Nprec = 1e-7; // Precision of e-foldings
const double dlogn = 0.1; // Width of bin in power spectrum
const double LL = 1.; // Box size L

double nsigmareset = 8.;
int aninum = 10; // dividing number for animation
const bool EoI_noise = false; // Adding the noise until EoN
int MeanNumber = 1; // Number of averageing in each
const bool sallslice = false; // save zeta for all zoom slice
int internumber = 0; // Number of zoom-in

// Outputs
const bool sfield = false; // Output field
const bool strajectory = false; // Output the trajectory // hybrid
const bool spower = true; // Output the power spectrum
bool sanimation = false; // Output the animation //hybrid
const bool sweight = false; // Output the weight
const bool scompaction = false; // Output the compaction

// Importance sampling
const double Nbias = 3.8; // Time of the bias
const double dNbias = 0.1; // Variance of the bias
const double bias = 0*sqrt(dNbias); // Amplitude of the bias

// Directory name of saved data, you can change after "make clean" in your terminal.
const std::string sdatadir = "data";
const std::string sourcedir = "../source";
const std::string noisefiledir = "noisedata/noiselist_";
const std::string biasfiledir = "biasdata/biaslist";
const std::string noisefilenamediv = "/part_";


const int NL = pow(NLnoise,3)/totalnoiseNo; // Box size L for each noisemap
double sqrt_dN = std::sqrt(dN);
constexpr int NLnoiseAll = NLnoise*NLnoise*NLnoise;
int totalstep = ceil(log((NLnoise/2-1)/sigma)/dN); // Total number of time step with noise
int firststep = ceil(log(nsigmareset/sigma)/dN);
int itpstep = ceil((log(nsigmareset/sigma)-log(nsigmareset/sigma/2.))/dN);
const double dx = LL/NLnoise; // Spacing of each lattice

#endif
