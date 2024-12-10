#ifndef INCLUDED_parameters_hpp_
#define INCLUDED_parameters_hpp_

#define MODEL 1 // 0 -> chaotic, 1 -> Starobinsky

// Parameters of STOLAS
const double sigma = 0.1; // ksigma = 2pi sigma exp(N) / L, nsigma = sigma exp(N)
const double dn = 1; // Thickness of nsigma sphere shell
const int NLnoise = pow(2,4); // Box size L
const int totalnoiseNo = pow(2,3); // The number of chunks
const int NL = pow(NLnoise,3)/totalnoiseNo; // Box size L for each noisemap
// const int divnumber = pow(2,5); // The number of chunks
const double dN = 0.01; // e-folds step
const double Nprec = 1e-7; // Precision of e-foldings
const double dlogn = 0.1; // Width of bin in power spectrum

// Importance sampling
const double Nbias = 3.8; // Time of the bias
const double dNbias = 0.1; // Variance of the bias
const double bias = 1*20*sqrt(dNbias); // Amplitude of the bias

// Outputs
const bool szeta = true; // Output the zeta
const bool spower = true; // Output the power spectrum
const bool scompaction = true; // Output the compaction
const bool sanimation = false; // Output the animation

// Directory name of saved data, you can change after "make clean" in your terminal.
const std::string sdatadir = "data";
const std::string sourcedir = "../source";
const std::string noisefilename = "noisedata/noisemap_";
const std::string biasfilename = "biasdata//biasmap.dat";
const std::string noisefiledir = "noisedata/noisetrs_";
const std::string noisefilenamediv = "/part_";
const std::string biasfilenamediv = "biasdata/biastrs/part_";

#endif
