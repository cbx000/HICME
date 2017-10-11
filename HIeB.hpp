#ifndef HIEBHPP

#define HIEBHPP

#include <string>

class HIeB
{
public:
  HIeB();
  double eBx, eBy, eBz; // 重离子碰撞中产生的磁场
  char flag; // 用来表明核运动方向的变量
  
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

  void CalVaccumEB(); // 计算真空磁场
  void CalQGPEB(); // 计算QGP中的磁场
  
private:
  double mx, my, mz, mt; // 计算磁场的时空坐标
  double meta, mtau; // 用eta和tau来表示z和t
  double mtau0; // QGP形成时间

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
