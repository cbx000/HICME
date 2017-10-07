#ifndef HIEBHPP

#define HIEBHPP

#include <string>

class HIeB
{
public:
  HIeB();
  double eBx, eBy, eBz; // 重离子碰撞中产生的磁场
  
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

  // 核参数
  std::string mNucleiType;
  double mR; // 核半径
  double md; // 核形状因子
  double mn0; // 归一化常数
  double mZ; // 核电荷数
  
  double ma; // 碰撞后快度分布参数

  double rhoFun_Ai(); // 核密度分布函数

};


#endif
