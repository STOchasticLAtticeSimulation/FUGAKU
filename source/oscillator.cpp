#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include "vec_op.hpp"
// #include "RK4.hpp" // RK4 を include

using namespace std;


template <class T>
void RK4(function<T(double, const T&)> dxdt, double &t, T &x, double dt);

// background x0, Fourier modes dx, その波数のリスト kk があるような系
template <class Tx0, class Tdx, class Tk>
void RK4(function<Tx0(double, const Tx0&)> dx0dt,
	 function<Tdx(double, const Tx0&, const Tdx&, const Tk&)> ddxdt,
	 double &t, Tx0 &x0, Tdx &dx, Tk &kk, double dt);


// t と x は更新するので参照渡し
template <class T>
void RK4(function<T(double, const T&)> dxdt, double &t, T &x, double dt)
{
  // ---------- RK4 の Butcher 係数 ------------------
  double a[4][4], b[4], c[4];

  a[0][0]=0;     a[0][1]=0;     a[0][2]=0;    a[0][3]=0;
  a[1][0]=1./2;  a[1][1]=0;     a[1][2]=0;    a[1][3]=0;
  a[2][0]=0;     a[2][1]=1./2;  a[2][2]=0;    a[2][3]=0;
  a[3][0]=0;     a[3][1]=0;     a[3][2]=1;    a[3][3]=0;
  
  b[0]=1./6;     b[1]=1./3;     b[2]=1./3;    b[3]=1./6;
  
  c[0]=0;        c[1]=1./2;     c[2]=1./2;    c[3]=1;
  // ------------------------------------------------


  T k[4]; // 4つの k_i
  T y; // k_i を求めるには y = x + dt sum_j^(i-1) a_ij k_j が必要


  for (int i=0; i<4; i++) {
    y = x;

    for (int j=0; j<i; j++) {
      y += dt * a[i][j] * k[j];
    }

    k[i] = dxdt(t+c[i]*dt, y);
  }

  
  t += dt; // t を更新

  for (int i=0; i<4; i++) {
    x += dt * b[i] * k[i]; // x を更新
  }
}


template <class Tx0, class Tdx, class Tk>
void RK4(function<Tx0(double, const Tx0&)> dx0dt,
	 function<Tdx(double, const Tx0&, const Tdx&, const Tk&)> ddxdt,
	 double &t, Tx0 &x0, Tdx &dx, Tk &kk, double dt)
{
  // ---------- RK4 の Butcher 係数 ------------------
  double a[4][4], b[4], c[4];

  a[0][0]=0;     a[0][1]=0;     a[0][2]=0;    a[0][3]=0;
  a[1][0]=1./2;  a[1][1]=0;     a[1][2]=0;    a[1][3]=0;
  a[2][0]=0;     a[2][1]=1./2;  a[2][2]=0;    a[2][3]=0;
  a[3][0]=0;     a[3][1]=0;     a[3][2]=1;    a[3][3]=0;
  
  b[0]=1./6;     b[1]=1./3;     b[2]=1./3;    b[3]=1./6;
  
  c[0]=0;        c[1]=1./2;     c[2]=1./2;    c[3]=1;
  // ------------------------------------------------


  Tx0 k0[4]; // x0 のための k_i
  Tdx dk[4]; // dx のための k_i
  Tx0 y0; // x0 のための y
  Tdx dy; // dx のための dy

  for (int i=0; i<4; i++) {
    y0 = x0;
    dy = dx;

    for (int j=0; j<i; j++) {
      y0 += dt * a[i][j] * k0[j];
      dy += dt * a[i][j] * dk[j];
    }

    k0[i] = dx0dt(t+c[i]*dt, y0);
    dk[i] = ddxdt(t+c[i]*dt, y0, dy, kk);
  }

  t += dt; // t を更新

  for (int i=0; i<4; i++) {
    x0 += dt * b[i] * k0[i]; // x0 を更新
    dx += dt * b[i] * dk[i]; // dx を更新
  }
}



vector<double> dxdt(double t, const vector<double> &x);

int main()
{
  double dt = 0.01;
  double t = 0;
  vector<double> x{1,0}; // 初期値
  double tf = 10; // 終了時刻

  ofstream ofs("oscillator.dat");
  ofs << std::setprecision(20);

  // t が tf に到達するまで RK4 を実行。
  while (t<tf) {
    // ターミナルとファイルに出力
    cout << t << ' ' << x[0] << ' ' << x[1] << endl;
    ofs << t << ' ' << x[0] << ' ' << x[1] << endl;
    
    RK4<vector<double>>(dxdt,t,x,dt);
    x[0] += 5.e-5;
  }
}


// 調和振動子 EoM
vector<double> dxdt(double t, const vector<double> &x)
{
  vector<double> dxdt(2);

  dxdt[0] = x[1];
  dxdt[1] = -x[0];

  return dxdt;
}