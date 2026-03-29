#ifndef INCLUDED_noise_bias_hpp_
#define INCLUDED_noise_bias_hpp_


// judge if point is in nsigma sphere shell
bool innsigma(int nx, int ny, int nz, int Num, double nsigma, double dn) {
  int nxt = (nx<=Num/2 ? nx : nx-Num);
  int nyt = (ny<=Num/2 ? ny : ny-Num);
  int nzt = (nz<=Num/2 ? nz : nz-Num);

  return std::abs(sqrt(nxt*nxt + nyt*nyt + nzt*nzt) - nsigma) <= dn/2.;
}

// judge real point
bool realpoint(int nx, int ny, int nz, int Num) {
  return (nx==0||nx==Num/2) && (ny==0||ny==Num/2) && (nz==0||nz==Num/2);
}

// judge independent complex point
bool complexpoint(int nx, int ny, int nz, int Num) {
  int nxt = (nx<=Num/2 ? nx : nx-Num);
  int nyt = (ny<=Num/2 ? ny : ny-Num);
  int nzt = (nz<=Num/2 ? nz : nz-Num);

  return (1<=nxt && nxt!=Num/2 && nyt!=Num/2 && nzt!=Num/2) ||
    (nxt==Num/2 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) || (nxt!=Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt!=Num/2) ||
    (nxt==0 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) ||
    (nxt==Num/2 && nyt==Num/2 && 1<=nzt && nzt!=Num/2) || (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt==Num/2) ||
    (nxt==0 && 1<=nyt && nyt!=Num/2 && nzt==0) ||
    (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==0) || (1<=nxt && nxt!=Num/2 && nyt==0 && nzt==Num/2) || (nxt==0 && nyt==Num/2 && 1<=nzt && nzt!=Num/2);
}

void dwlist_gen(double N, std::mt19937& engine) {
  int count = 0;
  double nsigma = sigma*exp(N);

  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] = 0.0;
    in[i][1] = 0.0;
  }
  
  LOOP{
    int idx = i*NLnoise*NLnoise + j*NLnoise + k;
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      if (realpoint(i,j,k,NLnoise)) {
	in[idx][0] = dist(engine);
	count++;
      } else if (complexpoint(i,j,k,NLnoise)) {
	in[idx][0] = dist(engine)/sqrt(2);
	in[idx][1] = dist(engine)/sqrt(2);
	count++;
      }
    }
  }

  // reflection
  int ip, jp, kp; // reflected index
  LOOP{
    int idx = i*NLnoise*NLnoise + j*NLnoise + k;
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      if (!(realpoint(i,j,k,NLnoise)||complexpoint(i,j,k,NLnoise))) {
        ip = (i==0 ? 0 : NLnoise-i);
        jp = (j==0 ? 0 : NLnoise-j);
        kp = (k==0 ? 0 : NLnoise-k);
	
	in[idx][0] = in[ip*NLnoise*NLnoise + jp*NLnoise + kp][0];
	in[idx][1] = -in[ip*NLnoise*NLnoise + jp*NLnoise + kp][1];
	count++;
      }
    }
  }

  if (count==0) {
    for (int i = 0; i < NLnoiseAll; i++) {
      dwlist[0][i] = out[i][0];
    }
    fftw_free(in);
    fftw_free(out);
    return;
  }

  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] /= sqrt(count);
    in[i][1] /= sqrt(count);
  }

  fftw_plan plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);

  for (int i = 0; i < NLnoiseAll; i++) {
    dwlist[0][i] = out[i][0];
  }

  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);
}

void biaslist1D(double N) {
  int count = 0;
  double nsigma = sigma*exp(N);
  
  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] = 0.0;
    in[i][1] = 0.0;
  }

  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      int idx = i*NLnoise*NLnoise + j*NLnoise + k;
      in[idx][0] = 1.0;
      in[idx][1] = 0.0;
      count++;
    }
  }

  if (count==0) {
    for (int i = 0; i < NLnoiseAll; i++) {
      biaslist[0][i] = out[i][0];
    }
    fftw_free(in);
    fftw_free(out);
    return;
  }
  for (int i = 0; i < NLnoiseAll; i++) {
    in[i][0] /= count;
  }

  fftw_plan plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);

  for (int i = 0; i < NLnoiseAll; i++) {
    biaslist[0][i] = out[i][0];
  }

  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);
}


#endif
