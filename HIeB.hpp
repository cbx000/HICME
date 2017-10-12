#ifndef HIEBHPP

#define HIEBHPP

#include <string>
#include <gsl/gsl_spline.h>

class HIeB
{
public:
  HIeB();
  double eBx, eBy, eBz; // 重离子碰撞中产生的磁场
  char flag; // 用来表明核运动方向的变量
  
  // 插值相关变量
  size_t N;
  double *ETA;
  double *EBY0;
  gsl_interp_accel * acc;
  gsl_spline *spline_steffen;

  // 设置场点的时空坐标
  void SetSpaceTime(double x, double y, double z, double t);
  // 设置场点的时空坐标(以另外一种坐标形式)
  void SetSpaceTime_tau(double tau, double eta);

  // 获取坐标
  double GetX() const;
  double GetY() const;
  double GetZ() const;
  double GetT() const;
  double GetEta() const;
  double GetTau() const;
  
  // 设置与获取质心系能量
  void SetSqrtS(double sqrtS);
  double GetSqrtS() const;
  // 设置与获取洛伦兹收缩因子gamma
  void SetGamma(double gamma);
  double GetGamma() const;
  // 设置与获取初始快度
  void SetY0(double Y0);
  double GetY0() const;

  // 设置与获取碰撞参数
  void SetB(double b);
  double GetB() const;

  // 设置与获取tau0
  void SetTau0(double tau0);
  double GetTau0() const;
  
  // 设置与获取核类型
  void SetNucleiType(std::string nuclei); // method to set nuclei parameters
  std::string GetNucleiType() const; // method to get nuclei name

  // 获取核参数
  double GetR() const; // 获取半径R
  double GetChargeZ() const; // 获取电荷量Z
  double GetN0() const; // 获取n0
  double GetD() const; // 获取d

  // 获取碰撞后快度分布参数a
  double GetA() const;
  
  // 设置与获取计算方法, ellipsoid为0, disklike为1
  void SetMethod(int method);
  int GetMethod() const;

  void CalVaccumEB(); // 计算磁场不考虑QGP响应
  void CaleBy00(); // 计算原点初始磁场
  void CaleBy0(size_t n);  // 计算沿z轴分布的初始磁场
  void CalOriginQGPeB(); // 计算原点磁场考虑QGP响应
  void CalQGPeB(); // 计算磁场考虑QGP响应
  
private:
  double mx, my, mz, mt; // 计算磁场的时空坐标
  double meta, mtau; // 用eta和tau来表示z和t
  double mtau0; // QGP形成时间

  int mIseBy00cal; // 是否已计算原点初始磁场, 0为否, 1为是
  int mIseBy0cal; // 是否已计算沿z轴的初始磁场, 0为否, 1为是
  double meBy00; // 原点初始磁场
  double mcs2, max2;

  double mb; // 碰撞参量
  double mGamma; // 洛伦兹收缩因子gamma
  double mSqrtS; // 质心系能量
  double mY0; // 初始快度

  int mMethod; // 计算方法, ellipsoid为0, disklike为1

  // 核参数
  std::string mNucleiType;
  double mR; // 核半径
  double md; // 核形状因子
  double mn0; // 归一化常数
  double mZ; // 核电荷数
  
  double ma; // 碰撞后快度分布参数

};

double rhoFun(double xp, double yp, double zp, char flag, double Y0, double b, double n0, double R, double d);
double f(double Y, double Y0, double a);
int eB_Part_Int(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
int eB_Spec_Int(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);

#endif
