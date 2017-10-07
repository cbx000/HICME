#include <cassert>
#include <iostream>
#include <string>
#include <cmath>
#include "HIeB.hpp"
#include "sq.h"

const double alpha_EM = 1.0 / 137.036; // fine-structure constant
const double m = 0.938272; // mass of proton
const double hbarc = 197.32696; // hbar * c
const double mpi = 140.0; // mass of pion
const double hbarcOverMpiSqure = (hbarc*hbarc)/(mpi*mpi);

HIeB::HIeB()
{
  // 设置默认坐标
  SetSpaceTime(0.0, 0.0, 0.0, 0.1);
  SetSqrtS(200.0); // 设置默认质心系能量
  
  ma = 0.5; // 设置快度分布参数
  SetB(8.0); // 设置默认碰撞参量
  
  SetNucleiType("Au"); // 设置默认核类型为Au
}

void HIeB::SetSpaceTime(double x, double y, double z, double t)
{
  mx = x;
  my = y;
  mz = z;
  mt = t;

  mtau = sqrt(t*t - z*z);
  meta = 0.5 * log( (t+z) / (t-z) );
}

void HIeB::SetSpaceTime_tau(double tau, double eta)
{
  mtau = tau;
  meta = eta;

  mx = 0.0;
  my = 0.0;
  mz = tau * sqrt(cosh(2*eta)*0.5 - 0.5);
  mt = tau * sqrt(cosh(2*eta)*0.5 + 0.5);
}

double HIeB::GetX() const
{
  return mx;
}

double HIeB::GetY() const
{
  return my;
}

double HIeB::GetZ() const
{
  return mz;
}

double HIeB::GetT() const
{
  return mt;
}

double HIeB::GetEta() const
{
  return meta;
}

double HIeB::GetTau() const
{
  return mtau;
}

void HIeB::SetSqrtS(double sqrtS)
{
  double E, p;
  
  assert (sqrtS >= 0.0);
  mSqrtS = sqrtS;

  // 计算对应的快度
  E = sqrtS / 2.0;
  p = sqrt(E*E - m*m);
  mY0 = 0.5 * log((E+p)/(E-p));

  // 计算对应的洛伦兹收缩因子
  mGamma = cosh(mY0);
}

double HIeB::GetSqrtS() const
{
  return mSqrtS;
}

void HIeB::SetGamma(double gamma)
{
  double p, V;
  assert (gamma >= 0.0);
  mGamma = gamma;

  // 计算对应的快度
  mY0 = acosh(gamma);

  // 计算对应的质心系能量
  V = tanh(mY0); // 速度(自然单位制)
  p = mGamma * m * V;
  mSqrtS = 2.0 * sqrt(m*m + p*p);
}

double HIeB::GetGamma() const
{
  return mGamma;
}

void HIeB::SetY0(double Y0)
{
  double p, V;
  assert (Y0 >= 0.0);
  mY0 = Y0;

  // 计算对应的洛伦兹收缩因子
  mGamma = cosh(mY0);

  // 计算对应的质心系能量
  V = tanh(mY0);
  p = mGamma * m * V;
  mSqrtS = 2.0 * sqrt(m*m + p*p);
}

double HIeB::GetY0() const
{
  return mY0;
}

void HIeB::SetB(double b)
{
  mb = b;
}

double HIeB::GetB() const
{
  return mb;
}

void HIeB::SetTau0(double tau0)
{
  mtau0 = tau0;
}

double HIeB::GetTau0() const
{
  return mtau0;
}


void HIeB::SetNucleiType(std::string nuclei)
{
  if (nuclei == "Au") {
    mNucleiType = "Au";
    mR = 6.38;
    md = 0.535;
    mn0 = 8.596268e-4;
    mZ = 79.0;
  } else if (nuclei == "Pb") {
    mNucleiType = "Pb";
    mR = 6.624;
    md = 0.549;
    mn0 = 7.69244e-4;
    mZ = 82.0;
  } else if (nuclei == "Cu"){
    mNucleiType = "Cu";
    mR = 4.214;
    md = 0.586;
    mn0 = 2.67894e-4;
    mZ = 29.0;
  } else {
    printf("Error: undefind nuclei type!");
  }
}

std::string HIeB::GetNucleiType() const
{
  return mNucleiType;
}

/*
double HIeB::rhoFun_Ai(char flag)
{
  double rho;
  double len;
  double newb = 0.0;

  // 判断是左核(+)还是右核(-)
  if (flag == '+')
    newb = -mb;
  else if (flag == '-')
    newb = mb;
  else
    std::cout << "Error: flag must be + or -!\n";

  len = sqrt(Sq(xp - newb/2.0) + Sq(yp) + Sq(mGamma*zp));
  rho = mGamma * mn0 / (1 + exp((len - mR)/md) ); // 此处乘以gamma是因为密度变大

  return rho;
}
*/

 /*
int eB_Part_Int(const int *ndim, const double xx[], const int *ncomp,
		      double ff[], void *userdata)
{
  double extra = 3.0; // extra is used to expand integral limits, because Woods-Saxon distribution is not precisely in the sphere of R
  double Imin[4]; // down limits
  double Imax[4]; // up limits

  Imin[0] = b/2.0 - mR;
  Imax[0] = -Imin[0];
  xp = Imin[0] + (Imax[0] - Imin[0]) * xx[0];

  Imin[1] = -sqrt(Sq(mR) - Sq(fabs(xp) + b/2.0));
  Imax[1] = -Imin[1];
  yp = Imin[1] + (Imax[1] - Imin[1]) * xx[1];
}
 */

 /*
int eB_Intgrand(const int *ndim, const double xx[], const int *ncomp,
		double ff[], void *userdata)
{
  int i;
  double sign, jacobian, extra, rho;
  double Imin[4];
  double Imax[4];

  extra = 3.0; // 加上extra是因为wood-saxon分布并不是完全在半径为R的球内

  // 根据被积区域类型和核标记确定积分上下限
  if (mRegionType == 'p') {
    Imin[0] = -(mR - b/2.0);
    Imax[0] = mR - b/2.0;
    Imin[1] = -sqrt(Sq(mR) - Sq(b/2.0));
    Imax[1] = sqrt(Sq(mR) - Sq(b/2.0));
    Imin[2] = -(mR + extra) / mGamma;
    Imax[2] = (mR + extra) / mGamma;
    Imin[3] = -mY0;
    Imax[3] = mY0;
  } else {
    // 判断核标记
    if (mFlag == '+') {
      Imin[0] = -(mR + b/2.0) - extra;
      Imax[0] = 0.0;
      Imin[1] = -mR - extra;
      Imax[1] = mR + extra;
      Imin[2] = -(mR+extra) / mGamma;
      Imax[2] = (mR+extra) / mGamma;
    } else {
      Imin[0] = 0.0;
      Imax[0] = (mR + b/2.0) + extra;
      Imin[1] = -mR - extra;
      Imax[1] = mR + extra;
      Imin[2] = -(mR + extra) / mGamma;
      Imax[2] = (mR + extra) / mGamma;
    }
  }

  // 变量变换 x -> min + (max - min) * x 将积分区间变为0-1
  xp = Imin[0] + (Imax[0] - Imin[0]) * xx[0];
  yp = Imin[1] + (Imax[1] - Imin[1]) * xx[1];
  zp = Imin[2] + (Imax[2] - Imin[2]) * xx[2];
  if (*ndim == 4)
    mY = Imin[3] + (Imax[3] - Imin[3]) * xx[3];

  if (mFlag == '+')
    sign = 1.0;
  else
    sign = -1.0;

  rho = rhoFun_Ai();
  
  if (mRegionType == 'p') {
    // 判断是否在被积区域内
    if ( (Sq(xp + b/2.0) + Sq(yp) <= Sq(mR)) &&
	 (Sq(xp - b/2.0) + Sq(yp) <= Sq(mR)) ) {
      pointEB_Wangqun(sign);
      ff[0] = mZ * rho * peBx;
      ff[1] = mZ * rho * peBy;
      ff[2] = mZ * rho * peBz;
    } else {
      ff[0] = 0.0;
      ff[1] = 0.0;
      ff[2] = 0.0;
    }
  } else if (mRegionType == 's') {
    if ( Sq(xp - sign*b/2.0) + Sq(yp) >= Sq(mR) ) {
      pointEB(sign);
      ff[0] = mZ * rho * peBx;
      ff[1] = mZ * rho * peBy;
      ff[2] = mZ * rho * peBz;
    } else {
      ff[0] = 0.0;
      ff[1] = 0.0;
      ff[2] = 0.0;
    }
  }

  jacobian = 1.0;
  for (i = 0; i < *ndim; i++) {
    jacobian = jacobian * (Imax[i] - Imin[i]);
  }

  ff[0] *= jacobian;
  ff[1] *= jacobian;
  ff[2] *= jacobian;

  return 0;
  
}
 */


