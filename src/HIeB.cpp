#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cuba.h>
#include "HIeB.hpp"
#include "sq.h"
using namespace std;

const double alpha_EM = 1.0 / 137.036; // 精细结构常数
const double m = 0.938272; // 质子质量
const double hbarc = 197.32696; // hbar * c
// const double mpi = 140.0; // π介子质量
// const double hbarcOverMpiSqure = (hbarc*hbarc)/(mpi*mpi);

HIeB::HIeB()
{
  // 设置默认坐标
  SetSpaceTime(0.0, 0.0, 0.0, 0.1);
  SetSqrtS(200.0); // 设置默认质心系能量
  SetB(8.0); // 设置默认碰撞参量
  SetTau0(0.1); // 设置默认tau0
  mIseBy00cal = 0; // 开始没有计算eBy00
  mIseBy0cal = 0; // 开始没有计算沿z轴的初始磁场
  mcs2 = 1.0/3.0;
  max2 = Sq(3.0);
  SetNucleiType("Au"); // 设置默认核类型为Au
  ma = 0.5; // 设置快度分布参数
  SetMethod(0); // 设置默认计算方法为ellipsoid
  SetLambda(0.2); // 设置屏蔽长度
  SetNpm(200.0);
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

  mIseBy0cal = 0;
  mIseBy00cal = 0;
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

  mIseBy0cal = 0;
  mIseBy00cal = 0;
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

  mIseBy0cal = 0;
  mIseBy00cal = 0;
}

double HIeB::GetY0() const
{
  return mY0;
}

void HIeB::SetB(double b)
{
  mb = b;
  mIseBy0cal = 0;
  mIseBy00cal = 0;
}

double HIeB::GetB() const
{
  return mb;
}

void HIeB::SetTau0(double tau0)
{
  mtau0 = tau0;
  mIseBy0cal = 0;
  mIseBy00cal = 0;
}

double HIeB::SetTau0byCen(double cen, double b)
{
  mb = b;// 之所以在函数参数中要加入b，是为了防止在通过中心度设置tau0时，忘记更新碰撞参数b

  double zeta = 2.0;
  double tau0_1 = zeta * mR * exp(-mY0);

  double lambda = 0.3; // 注意：这个参数与mlambda的含义完全不一样
  double A_0 = 197.0;
  double sqrts_0 = 130.0;
  double b_0 = 1.611 * pow(cen, 0.4817) - 0.2972;
  double Qs2_0 = -0.0004779 * pow(b_0, 3.0) - 0.004795 * Sq(b_0) - 0.005726*b_0+2.051;
  double Qs2 = Qs2_0 * pow(mA/A_0,2.0/(6.0+3.0*lambda)) * pow(mSqrtS/sqrts_0,2.0*lambda/(2.0+lambda));
  double tau0_2 = 1.0/sqrt(Qs2)*0.1975; // 0.1975 是单位转换因子 
  if (tau0_2 > tau0_1) {
    mtau0 = tau0_2;
  } else {
    mtau0 = tau0_1;
  }
  return mtau0;
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
    mA = 197.0;
  } else if (nuclei == "Pb") {
    mNucleiType = "Pb";
    mR = 6.624;
    md = 0.549;
    mn0 = 7.69244e-4;
    mZ = 82.0;
    mA = 207.0;
  } else if (nuclei == "Cu"){
    mNucleiType = "Cu";
    mR = 4.214;
    md = 0.586;
    mn0 = 2.67894e-4;
    mZ = 29.0;
    mA = 63.0;
  } else {
    printf("Error: undefind nuclei type!");
  }

  mIseBy0cal = 0;
  mIseBy00cal = 0;
}

std::string HIeB::GetNucleiType() const
{
  return mNucleiType;
}

double HIeB::GetR() const
{
  return mR;
}

double HIeB::GetChargeZ() const
{
  return mZ;
}

double HIeB::GetN0() const
{
  return mn0;
}

double HIeB::GetD() const
{
  return md;
}

double HIeB::GetA() const
{
  return ma;
}

void HIeB::SetMethod(int method)
{
  mMethod = method;

  mIseBy0cal = 0;
  mIseBy00cal = 0;
}

int HIeB::GetMethod() const
{
  return mMethod;
}


void HIeB::CalVaccumEB()
{
  double totalerror;
  int neval, fail;
  double interror, prob;
  double eBp_plus, eBp_minus, eBs_plus, eBs_minus;

  int nvec = 1;
  double epsrel = 0.01;
  double epsabs = 0.5;
  int flags = 0 | 8;
  int seed = 0;
  int smineval = 1.5e4;
  int smaxeval = 4e6;
  int pmineval = 2e4;
  int pmaxeval = 1e7;
  // int mineval = 1e5;
  // int maxeval = 1e7;
  int nstart = 1000;
  int nincrease = 500;
  int nbatch = 1000;
  int gridno = 0;
  char *statefile = NULL;
  void *spin = NULL;
  
  totalerror = 0.0;
  // 对于“参与者”
  flag = '+'; // 正向
  Vegas(4, 1, eB_Part_Int, (void *)this, nvec, epsrel, epsabs, flags, seed,
        pmineval, pmaxeval, nstart, nincrease, nbatch, gridno,
        statefile, spin, &neval, &fail, &eBp_plus, &interror, &prob);
  totalerror += interror;
  flag = '-'; // 负向
  Vegas(4, 1, eB_Part_Int, (void *)this, nvec, epsrel, epsabs, flags, seed,
        pmineval, pmaxeval, nstart, nincrease, nbatch, gridno,
        statefile, spin, &neval, &fail, &eBp_minus, &interror, &prob);
  totalerror += interror;
  
  // 对于“旁观者”
  flag = '+'; // 正向
  Vegas(3, 1, eB_Spec_Int, (void *)this, nvec, epsrel, epsabs, flags, seed,
        smineval, smaxeval, nstart, nincrease, nbatch, gridno,
        statefile, spin, &neval, &fail, &eBs_plus, &interror, &prob);
  totalerror += interror;
  flag = '-'; // 负向
  Vegas(3, 1, eB_Spec_Int, (void *)this, nvec, epsrel, epsabs, flags, seed,
        smineval, smaxeval, nstart, nincrease, nbatch, gridno,
        statefile, spin, &neval, &fail, &eBs_minus, &interror, &prob);
  totalerror += interror;
  eBy = eBp_plus + eBp_minus + eBs_plus + eBs_minus;
  // printf(" eBp_plus = %-8g\n eBp_minus = %-8g\n eBs_plus = %-8g\n eBs_minus = %-8g\n", eBp_plus, eBp_minus, eBs_plus, eBs_minus);

}


void HIeB::CaleBy00()
{
  double x, y, z, t;
  x = mx;
  y = my;
  z = mz;
  t = mt;
  SetSpaceTime(0.0, 0.0, 0.0, mtau0);
  CalVaccumEB();
  meBy00 = eBy;
  mIseBy00cal = 1;

  SetSpaceTime(x, y, z, t);
}

double HIeB::GeteBy00()
{
  if (mIseBy00cal != 1) {
    CaleBy00();
  }

  return meBy00;
}

void HIeB::CaleBy0(int n)
{
  double minEta, maxEta;
  double x = mx, y = my, z = mz, t = mt;
  minEta = -mY0-0.5;
  maxEta = mY0+0.5;

  N = n+1;

  ETA = (double *)malloc((N) * sizeof(double));
  EBY0 = (double *)malloc((N) * sizeof(double));

  string filename = "data/" + mNucleiType + "_" + to_string(mSqrtS) 
    + "_"+ to_string(mb) + "_"+ to_string(mMethod) +"_"+ to_string(N) +".dat";

  if (ifstream(filename)) { // 如果文件存在，则从文件中读取数据
    ifstream myfile (filename);
    if (myfile.is_open()) {
      for (int i = 0; i <= n; i++) {
        myfile >> ETA[i];
        myfile >> EBY0[i];
      }
      myfile.close();
    }
  } else { // 如果文件不存在，则计算数据，并保存到文件中
    ofstream myfile (filename);
    if (myfile.is_open()) {
      for (int i = 0; i <= n; i++) {
        printf("\r 正在计算初始磁场: %3.0lf %%", (double)i/n*100.0);
        ETA[i] = minEta + (double)i*(maxEta - minEta)/(double)N;
        SetSpaceTime_tau(mtau0, ETA[i]);
        CalVaccumEB();
        EBY0[i] = eBy;
        myfile << ETA[i] << " " << EBY0[i] << endl;
      }
      myfile.close();
    }  else {
      cout << "无法打开文件" << endl;
      return;
    }
  }

  spline_steffen = gsl_spline_alloc(gsl_interp_steffen, N);
  gsl_spline_init(spline_steffen, ETA, EBY0, N);

  SetSpaceTime(x, y, z, t);
  mIseBy0cal = 1;
}

void HIeB::CalOriginQGPeB()
{
  if (mIseBy00cal != 1) {
    CaleBy00();
  }
  
  eBy = (mtau0 / mtau) * exp(- mcs2/(2.0*max2)*(Sq(mtau) - Sq(mtau0))) * meBy00;
    
}

void HIeB::CalQGPeB()
{
  if (mIseBy0cal != 1) {
    HIeB::CaleBy0(100);
  }

  double cosheta = cosh(meta);
  eBy = mtau0/mtau*exp(-mcs2/(2.0*max2)*(Sq(mtau) - Sq(mtau0))*Sq(cosheta))*gsl_spline_eval(spline_steffen, meta, acc);
}

double HIeB::xifun(double xp, double yp, char sign)
{
  double y_plus, y_minus;
  y_plus = sqrt(Sq(mR) - Sq(fabs(xp) + mb/2.0));
  y_minus = - y_plus;
  if (sign == '+') {
    return exp(-fabs(y_plus - yp)/mlambda);
  } else {
    return exp(-fabs(y_minus - yp)/mlambda);
  }
}

void HIeB::SetLambda(double lambda)
{
  mlambda = lambda * mR;

  mIseBy0cal = 0;
  mIseBy00cal = 0;
}

double HIeB::GetLambda()
{
  return mlambda;
}

void HIeB::SetNpm(double npm)
{
  mNm = npm;
  mNp = npm;
}

double HIeB::GetNp()
{
  return mNp;
}

double HIeB::GetNm()
{
  return mNm;
}

void HIeB::cmefun()
{
  int nvec = 1, neval, flags = 0 | 8, seed = 0, fail;
  int mineval = 1e5, maxeval = 1e7, nstart = 1000, nincrease = 500, nbatch = 1000, gridno = 0;
  char *statefile = NULL;
  void *spin = NULL;
  double epsrel = 0.01, epsabs = 0.5, interror, prob;
  
  if (mIseBy00cal != 1) {
    CaleBy00();
  }

  Vegas(3,1,delta_pp_Int, (void *)this, nvec, epsrel, epsabs, flags, seed, mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile, spin, &neval, &fail, &delta_pp, &interror, &prob);
  Vegas(3,1,delta_pm_Int, (void *)this, nvec, epsrel, epsabs, flags, seed, mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile, spin, &neval, &fail, &delta_pm, &interror, &prob);
  printf("delta_pp = %g delta_pm = %g\n", delta_pp, delta_pm);

  app = 1.0/Sq(mNp)*Sq(M_PI)/16.0*(delta_pp);
  apm = 1.0/(mNp*mNm)*Sq(M_PI)/16.0*(delta_pm);
  printf("a_pp = %g a_pm = %g\n", app, apm);
}

double rhoFun(double xp, double yp, double zp, char flag, double Y0, double b, double n0, double R, double d) {
  // xp, yp, zp: 源点坐标
  // 返回: (xp, yp, zp) 处的数密度

  double rho;
  double len;
  double gamma;
  double sign = 1.0;
  gamma = cosh(Y0);
  
  if (flag == '+') {
    sign = 1.0; // 左核
  } else if (flag == '-') {
    sign = -1.0; // 右核
  } else {
    printf("[rhoFun]error: nucleus type(flag) should be '+' or '-'\n");
  }

  len = sqrt(Sq(xp + sign*b/2.0) + Sq(yp) + Sq(gamma*zp));
  rho = gamma * n0 / ( 1.0 + exp((len-R)/d) ); // 此处乘以gamma是因为密度增加了
  
  return rho;
}

double f(double Y, double Y0, double a) {
  // 参与者的快度分布函数
  return (a*exp(a*Y)) / (2*sinh(a*Y0));
}

int eB_Part_Int(const int *ndim, const double xx[], const int *ncomp,
			   double ff[], void *userdata) {
  // ndim是指向积分维度的指针
  // xx是指向积分变量数组的指针
  // ncomp是指向积分分量数的指针
  // ff是指向积分值数组的指针
  // userdata是指向传递参数数据的指针
  HIeB *ud = (HIeB *) userdata;
  double x = ud->GetX();
  double y = ud->GetY();
  double z = ud->GetZ();
  double t = ud->GetT();
  double Y0 = ud->GetY0();
  double b = ud->GetB();
  double R = ud->GetR();
  double Z = ud->GetChargeZ();
  char flag = ud->flag;
  int method = ud->GetMethod();
  double a = ud->GetA();
  double n0 = ud->GetN0();
  double d = ud->GetD();
  
  static double sign;
  static double jacobian;
  static double denominator;

  static double xp; // 积分变量: x'
  static double yp; // 积分变量: y'
  static double zp; // 积分变量: z'
  static double Y;  // 积分变量: Y

  static double Imin[4]; // 积分下限数组
  static double Imax[4]; // 积分上限数组
  static double gamma; // 洛伦兹收缩因子
  static double extra;

  gamma = cosh(Y0);
  extra = 3; // extra 是用来略微扩大积分限的，因为 Wood-Saxon 分布不是完全分布在半径为R的球内

  Imin[0] = b/2.0 - R; 
  Imax[0] = -Imin[0];
  // x'的上下限: b/2-R <= x' <= -b/2+R
  xp = Imin[0] + (Imax[0] - Imin[0]) * xx[0]; // 对xx[0]进行放缩

  Imin[1] = - sqrt(Sq(R) - Sq(fabs(xp) + b/2.0)); 
  Imax[1] = -Imin[1];
  // y'的上下限: -\sqrt{R^2 - (|x| + b/2)^2} <= y' <= \sqrt{R^2 - (|x'| + b/2)^2}
  yp = Imin[1] + (Imax[1] - Imin[1]) * xx[1]; // 对xx[1]进行放缩

  Imin[2] = -(R+extra)/(gamma); 
  Imax[2] = -Imin[2];
  // z'的上下限: -R/cosh Y_0 <= z' <=  R/cosh Y_0, 注意gamma = cosh(Y_0)
  zp = Imin[2] + (Imax[2] - Imin[2]) * xx[2]; // 对xx[2]进行放缩

  Imin[3] = -Y0;
  Imax[3] = Y0;
  // Y的上下限: -Y_0 <= Y <= Y_0
  Y = Imin[3] + (Imax[3] - Imin[3]) * xx[3]; // 对xx[3]进行放缩

  if (flag == '+') 
     sign = 1.0; // 沿z轴正向
  else
     sign = -1.0; // 沿z轴负向

  double temp = 0.0;
  if (method == 0) {
    temp = (zp/tanh(Y0) + sign * t )*sinh(Y) - z * cosh(Y);
  } else if (method == 1) {
    temp = ( sign * t )*sinh(Y) - z * cosh(Y);
  } else {
    printf("Error: method must be 0(ellipsoid) or 1(disklike)");
  }
  // temp是分母中的第二项
  
  denominator = pow(Sq(xp - x) + Sq(yp - y) +
                    Sq(temp),1.5);
  ff[0] = sign * Sq(hbarc) * Z * alpha_EM * f(Y, Y0, a) * sinh(Y) * rhoFun(xp, yp, zp, flag, Y0, b, n0, R, d) * (x - xp) / denominator;
  // 计算被积函数值, 注意Sq(hbarc)是用来转换单位的

  jacobian = 1.0;
  for (int i = 0; i < *ndim; i++) {
    jacobian = jacobian * (Imax[i] - Imin[i]);
  }
  ff[0] = jacobian * ff[0];
  // 由于积分变量被放缩了，所以结果必须乘以Jacobian
  
  return 0;
}

int eB_Spec_Int(const int *ndim, const double xx[], const int *ncomp,
			   double ff[], void *userdata) {
  // ndim是指向积分维度的指针
  // xx是指向积分变量数组的指针
  // ncomp是指向积分分量数的指针
  // ff是指向积分值数组的指针
  // userdata是指向传递参数数据的指针
  
  HIeB *ud = (HIeB *) userdata;
  double x = ud->GetX();
  double y = ud->GetY();
  double z = ud->GetZ();
  double t = ud->GetT();
  double Y0 = ud->GetY0();
  double b = ud->GetB();
  double R = ud->GetR();
  double Z = ud->GetChargeZ();
  char flag = ud->flag;
  int method = ud->GetMethod();
  double n0 = ud->GetN0();
  double d = ud->GetD();
  
  static double sign;
  static double jacobian;
  static double denominator;

  static double phip; // 积分变量phi'
  static double xpperp; // 积分变量x'_\perp
  static double zp; // 积分变量z'
  static double xp; 
  static double yp; 

  static double Imin[3]; // 积分下限数组
  static double Imax[3]; // 积分上限数组
  static double gamma; // 洛伦兹收缩因子
  static double extra;

  gamma = cosh(Y0);
  extra = 3; // extra 是用来略微扩大积分限的，因为 Wood-Saxon 分布不是完全分布在半径为R的球内

  if (flag == '+') { // 沿z轴正向
     sign = 1.0; 
     Imin[0] = M_PI/2.0;
     Imax[0] = 3.0*M_PI/2.0;
     // 对正向运动的来说，phi'的上下限: pi/2 <= phi' <= 3pi/2
     phip = Imin[0] + (Imax[0] - Imin[0]) * xx[0]; // 放缩xx[0]
  } else { // 沿z轴负向
     sign = -1.0; 
     Imin[0] = -M_PI/2.0;
     Imax[0] = M_PI/2.0;
     // 对负向运动的来说，phi'的上下限: -pi/2 <= phi' <= pi/2
     phip = Imin[0] + (Imax[0] - Imin[0]) * xx[0]; // 放缩xx[0]
  }

  Imin[1] = -b/2.0 * fabs(cos(phip)) + sqrt(Sq(R) - Sq(b)/4.0*Sq(sin(phip)));
  Imax[1] = b/2.0 * fabs(cos(phip)) + sqrt(Sq(R) - Sq(b)/4.0*Sq(sin(phip)));
  // x'_\perp的上下限: $-\frac{b}{2}cos(\phi') + \sqrt{R^2 - \frac{b^2}{4} \sin^2(\phi')} \leq x'_\perp \leq \frac{b}{2}cos(\phi') + \sqrt{R^2 - \frac{b^2}{4} \sin^2(\phi')}$
  xpperp = Imin[1] + (Imax[1] - Imin[1]) * xx[1]; // 放缩xx[1]

  Imin[2] = -(R+extra)/(gamma); 
  Imax[2] = -Imin[2];
  // z'的上下限: $-R/\cosh Y_0 \leq z' \leq  R/\cosh Y_0$, note that $\gamma = \cosh(Y_0)$ 
  zp = Imin[2] + (Imax[2] - Imin[2]) * xx[2]; // 放缩xx[2]

  xp = xpperp * cos(phip);
  yp = xpperp * sin(phip);

  double temp = 0.0;
  if (method == 0) {
    temp = (zp/tanh(Y0) + sign * t )*sinh(Y0) - z * cosh(Y0);
  } else if (method == 1) {
    temp = ( sign * t )*sinh(Y0) - z * cosh(Y0);
  } else {
    printf("Error: method must be 0(ellipsoid) or 1(disklike)");
  }
  // temp是分母的第二项
  
  denominator = pow(Sq(xp - x) + Sq(yp - y) +
                    Sq(temp),1.5);
  ff[0] = sign * Sq(hbarc) * Z * alpha_EM * sinh(Y0) * xpperp * rhoFun(xp, yp, zp, flag, Y0, b, n0, R, d) * (x - xp) / denominator;
  // 计算被积函数值, 注意Sq(hbarc)是用来转换单位的
  

  jacobian = 1.0;
  for (int i = 0; i < *ndim; i++) {
    jacobian = jacobian * (Imax[i] - Imin[i]);
  }
  ff[0] = jacobian * ff[0];
  // 由于积分变量被放缩了，所以结果必须乘以Jacobian
  return 0;
}

int delta_pp_Int(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
  HIeB *ud = (HIeB *) userdata;
  static double b, R, t0;
  b = ud->GetB();
  R = ud->GetR();
  t0 = ud->GetTau0();
  
  static double jacobian;
  static double xp;
  static double yp;
  static double tau;

  static double Imin[3];
  static double Imax[3];
  static double result;
  static double kappa = 1.0, alpha_s = 1.0, sumqf22 = 25.0/81.0;
  static double xi_plus, xi_minus;

  static double tempeB;
  
  Imin[0] = b/2.0 - R;
  Imax[0] = -Imin[0];
  xp = Imin[0] + (Imax[0] - Imin[0]) * xx[0];
  Imin[1] = -sqrt(Sq(R) - Sq(fabs(xp) + b/2.0));
  Imax[1] = -Imin[1];
  yp = Imin[1] + (Imax[1] - Imin[1]) * xx[1];
  Imin[2] = t0; // 2.0 * R * exp(-Y0);
  Imax[2] = 15.0;
  tau = Imin[2] + (Imax[2] - Imin[2]) * xx[2];

  xi_plus = ud->xifun(xp, yp, '+');
  xi_minus = ud->xifun(xp, yp, '-');

  ud->SetSpaceTime_tau(tau, 0.0);
  ud->CalOriginQGPeB();
  tempeB = ud->eBy/Sq(hbarc); // 注意这里要转换单位
  result = 2.0 * kappa * alpha_s * sumqf22 * (Sq(xi_plus) + Sq(xi_minus)) * tau * Sq(tempeB);

  jacobian = 1.0;
  for (int i = 0; i < *ndim; i++) {
    jacobian = jacobian * (Imax[i] - Imin[i]);
  }
  ff[0] = jacobian * result;
  
  return 0;
}

int delta_pm_Int(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
  HIeB *ud = (HIeB *)userdata;
  static double b, R, t0;
  b = ud->GetB();
  R = ud->GetR();
  t0 = ud->GetTau0();
  
  static double jacobian;
  static double xp;
  static double yp;
  static double tau;

  static double Imin[3];
  static double Imax[3];
  static double result;
  static double kappa = 1.0, alpha_s = 1.0, sumqf22 = 25.0/81.0;
  static double xi_plus, xi_minus;
  
  static double tempeB;

  Imin[0] = b/2.0 - R;
  Imax[0] = -Imin[0];
  xp = Imin[0] + (Imax[0] - Imin[0]) * xx[0];
  Imin[1] = -sqrt(Sq(R) - Sq(fabs(xp) + b/2.0));
  Imax[1] = -Imin[1];
  yp = Imin[1] + (Imax[1] - Imin[1]) * xx[1];
  Imin[2] = t0; // 2.0 * R * exp(-Y0);
  Imax[2] = 15.0;
  tau = Imin[2] + (Imax[2] - Imin[2]) * xx[2];

  xi_plus = ud->xifun(xp, yp, '+');
  xi_minus = ud->xifun(xp, yp, '-');

  ud->SetSpaceTime_tau(tau, 0.0);
  ud->CalOriginQGPeB();
  tempeB = ud->eBy/Sq(hbarc); // 注意这里要转换单位
  result = -4.0 * kappa * alpha_s * sumqf22 * xi_plus * xi_minus * tau * Sq(tempeB);

  jacobian = 1.0;
  for (int i = 0; i < *ndim; i++) {
    jacobian = jacobian * (Imax[i] - Imin[i]);
  }
  ff[0] = jacobian * result;
  
  return 0;
}
