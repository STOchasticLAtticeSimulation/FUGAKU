#ifndef INCLUDED_output_hpp_
#define INCLUDED_output_hpp_


void OpenFiles(int NoisefiledirNo, int Interpolatingnumber){
  std::string NLfilename = std::to_string(NLnoise) + std::string("_") + std::to_string(NFIELDS) + std::string("_") + std::to_string(NoisefiledirNo);
  std::string InterFileName = NLfilename + std::string("_") + std::to_string(Interpolatingnumber);

  // Nfile.open(Nfileprefix + InterFileName + std::string(".dat"));
  Nfile.open(Nfileprefix + InterFileName + std::string(".bin"), std::ios::binary);
  Nfile << std::setprecision(10);
  Nfilefail = Nfile.fail();

  if (sfield) {
    fieldfile.open(fieldfileprefix + InterFileName + std::string(".dat"));
    fieldfile << std::setprecision(10);
  }
  
  if (strajectory) {
    trajectoryfile.open(trajectoryfileprefix + InterFileName + std::string(".dat"));
    trajectoryfile << std::setprecision(10);
  }

  if (sweight) {
    logwfile.open(logwfileprefix + InterFileName + std::string(".dat"));
  }
  
}

void save_N1N2(int NoisefiledirNo){
  std::ofstream N1N2fileA;
  std::string fFileName = std::to_string(NLnoise) + std::string("_") + std::to_string(NFIELDS) + std::string("_") + std::to_string(NoisefiledirNo) + std::string("_") + std::to_string(0);
  N1N2fileA.open(sdatadir + "/" + model + "/N1N2_" + fFileName + std::string(".dat"));
  N1N2fileA << std::setprecision(10);

  double N1av=0., N2av=0.;
  for (int i=0; i<NLnoiseAll; i++){
    N1av += N1list[i];
    N2av += N2list[i];
  }  
  N1N2fileA << N1av/NLnoiseAll << ' ' << N2av/NLnoiseAll << std::endl;
  std::cout << "Export N1 and N2" << std::endl;
  N1N2fileA.close();
}

void save_zeta(){
  // for (int i=0; i<NLnoiseAll; i++) Nfile << i << ' ' << Ndata[i] << std::endl;
  Nfile.write(reinterpret_cast<const char*>(&Ndata), sizeof(double) * NLnoiseAll);
  std::cout << "Export delta N map" << std::endl;
}

void save_field(){
  for (int i=0; i<NLnoiseAll; i++) {
    #if MODEL==3
      NFLOOP{
        fieldfile << phievol[i][2*nf] << ' ';
      }
      fieldfile << std::endl;
    #else
      fieldfile << phievol[i][0] << std::endl;
    #endif
  }
}

void save_trajectory(state_type PHI, double Ntime){
  trajectoryfile << Ntime << ' ' << PHI[0];
  #if MODEL==3
    double PSISUM = 0;
    NFLOOP PSISUM += pw2(PHI[2*nf]);
    trajectoryfile << ' ' << PSISUM << ' ' << EoI(PHI)-2 << ' ' << hubble(PHI);
  #endif
  trajectoryfile << std::endl;
}

void animation(std::array<state_type,NLnoiseAll>& phievol, int NoisefiledirNo, int Ntime){
  std::string fFileName = std::to_string(NLnoise) + std::string("_") + std::to_string(NFIELDS) + std::string("_") + std::to_string(NoisefiledirNo) + std::string("_") + std::to_string(Ntime);
  // fieldfileA.open(animationfileprefix + fFileName + std::string(".dat"));
  fieldfileA.open(animationfileprefix + fFileName + std::string(".bin"), std::ios::binary);
  fieldfileA << std::setprecision(10);

  static std::array<double, 2 * NLnoiseAll> buffer;

  for (int i=0; i<NLnoiseAll; i++) {
      buffer[2*i] = phievol[i][1];
      buffer[2*i + 1] = hubble(phievol[i]);
  }
  fieldfileA.write(reinterpret_cast<const char*>(buffer.data()), sizeof(double) * buffer.size());

  // for (int i=0; i<NLnoiseAll; i++){
  //   fieldfileA << phievol[i][1] << ' ' << hubble(phievol[i]) << std::endl;

  //   #if MODEL==3
  //     double PSISUM = 0;
  //     NFLOOP PSISUM += pw2(phievol[i][2*nf]);
  //     fieldfileA << sqrt(PSISUM) << std::endl;

  //   #endif
  // }
  fieldfileA.close();
}


// Calculate power spectrum
void spectrum(std::array<double,NLnoiseAll>& Ndata, int noisefiledirNo) {  
  powsfile.open(powsfileprefix + std::string(".dat"), std::ios::app);
  powsfile << std::setprecision(10);

  fft_1D_real(Ndata);

  LOOP{
    int nxt, nyt, nzt; // shifted index
    nxt = (i<=NLnoise/2 ? i : i-NLnoise);
    nyt = (j<=NLnoise/2 ? j : j-NLnoise);
    nzt = (k<=NLnoise/2 ? k : k-NLnoise);
    int idx = i*NLnoise*NLnoise + j*NLnoise + k;
    
    double rk=nxt*nxt+nyt*nyt+nzt*nzt;

    double LogNk = log(sqrt(rk));
    double calPk = norm(bkspectrum[idx])/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise;
    for (size_t ii = 0; ii < imax; ii++) {
      if ((dlogn*ii-LogNk)<dlogn/2. && -dlogn/2.<=(dlogn*ii-LogNk) && rk!=0) {
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
  powsfile.close();
  std::cout << "ExportPowerSpectrum" << std::endl;
}


// calculation of weight
void weight(int seed) {
  double logw = 0.;
  for (size_t n=0; n<weightlist.size(); n++) {
    double N = n*dN;
    double Bias = bias/dNbias/sqrt(2*M_PI)*exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
    logw -= Bias*weightlist[n]*sqrt(dN) + (Bias*Bias*dN)/2;
  }
  logwfile << seed << ' ' << logw << std::endl;
}


// Calculate compaction function
void compaction(std::array<double,NLnoiseAll>& Ndata, int noisefiledirNo) {
  prbfile.open(prbfileprefix + std::string(".dat"), std::ios::app);
  cmpfile.open(cmpfileprefix + std::to_string(NLnoise) + std::string("_") + std::to_string(noisefiledirNo) + std::string(".dat"));
  prbfile << std::setprecision(10);
  cmpfile << std::setprecision(10);
  
  double Naverage = 0;
  int dr = 1;
  for (size_t n = 0; n < Ndata.size(); n++) {
    Naverage += Ndata[n];
  }
  Naverage /= NLnoise*NLnoise*NLnoise;

  // zeta map
  for (size_t n = 0; n < Ndata.size(); n++) {
    Ndata[n] -= Naverage;
  }

  // Find max value
  int maxNpoint = std::distance(Ndata.begin(), std::max_element(Ndata.begin(), Ndata.end()));
  // int xmax = maxNpoint/NLnoise/NLnoise, ymax = (maxNpoint%(NLnoise*NLnoise))/NLnoise, zmax = maxNpoint%NLnoise;
  int xmax = 0, ymax = 0, zmax = 0;

  // radial profile
  for (size_t i=0; i<NLnoise*NLnoise*NLnoise; i++) {
    int nx=i/NLnoise/NLnoise ,ny=(i%(NLnoise*NLnoise))/NLnoise, nz=i%NLnoise;
    
    // centering
    int nxt, nyt, nzt; // shifted index
    nxt = (nx<=NLnoise/2+xmax ? nx-xmax : nx-xmax-NLnoise);
    nyt = (ny<=NLnoise/2+ymax ? ny-ymax : ny-ymax-NLnoise);
    nzt = (nz<=NLnoise/2+zmax ? nz-zmax : nz-zmax-NLnoise);

    for (size_t ri=0; ri<NLnoise/2; ri++) {
      double norm = std::abs(sqrt(nxt*nxt+nyt*nyt+nzt*nzt)-ri);
      if (norm<=dr/2.) {
        zetar[0][ri]++;
        zetar[1][ri]+=Ndata[i];
        break;
      }
    }
  }
  for (size_t ri=0; ri<NLnoise/2; ri++) {
    zetar[1][ri] /= zetar[0][ri]; // average
  }

  // derivative zeta
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

#endif