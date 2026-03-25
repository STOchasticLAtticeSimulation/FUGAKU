#ifndef INCLUDED_laplacian_hpp_
#define INCLUDED_laplacian_hpp_

inline int INCREMENT(int i) {
  return (i == NLnoise - 1) ? 0 : i + 1;
}

inline int DECREMENT(int i) {
  return (i == 0) ? NLnoise - 1 : i - 1;
}

std::vector<double> computeLaplacian(const std::vector<double> &x)
{
  std::vector<double> laplacian(NLnoise, 0.0);

  LOOP
  {
    int idx = index(i, j, k);
    int ip1 = INCREMENT(i);
    int im1 = DECREMENT(i);
    int jp1 = INCREMENT(j);
    int jm1 = DECREMENT(j);
    int kp1 = INCREMENT(k);
    int km1 = DECREMENT(k);

    laplacian[idx] = x[index(ip1, j, k)] + x[index(im1, j, k)] +
                     x[index(i, jp1, k)] + x[index(i, jm1, k)] +
                     x[index(i, j, kp1)] + x[index(i, j, km1)] -
                     6. * x[idx];
  }
  double norm = 1. / dx / dx;

  return norm * laplacian;
}

#endif
