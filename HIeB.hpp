#ifndef HIEBHPP

#define HIEBHPP

#include <cmath>
#include <string>

class HIeB
{
public:
  HIeB();
  double eBx, eBy, eBz; // 重离子碰撞中产生的磁场
  double x, y, z, t; // 计算磁场的时空坐标
  double b; // 碰撞参量
  double sigma, sigmaChi; // 电导率和手征磁导率(目前用不到)

  // 设置与获取质心系能量
  void SetSqrtS(double sqrtS);
  double GetSqrtS() const;
  // 设置与获取洛伦兹收缩因子gamma
  void SetGamma(double gamma);
  double GetGamma() const;
  // 设置与获取初始快度
  void SetY0(double Y0);
  double GetY0() const;
  // 设置与获取核类型
  void SetNucleiType(std::string nuclei); // method to set nuclei parameters
  std::string GetNucleiType() const; // method to get nuclei name

  void CalVaccumEB(int comp); // 计算真空磁场
  void CalQGPEB(int comp); // 计算QGP中的磁场
  
private:
  double mGamma; // 洛伦兹收缩因子gamma
  double mSqrtS; // 质心系能量
  double mY0; // 初始快度

  // 核参数
  std::string mNucleiType;
  double mR; // 核半径
  double md; // 核形状因子
  double mn0; // 归一化常数
  double mZ; // 核电荷数
  
  double ma; // 碰撞后快度分布参数


  char mRegionType, mFlag;

  double rhoFun_Ai();

  int eB_Part_Int(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
  int eB_Spec_Int(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
  int eB_Intgrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
  
};


#endif
