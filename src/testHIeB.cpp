#include <iostream>
#include <iomanip>
#include <cmath>
#include "HIeB.hpp"
using namespace std;

int main(int argc, char *argv[])
{
  double sqrtS, gamma, Y0, tolerance = 1e-3, t;
  HIeB myeB;

  setprecision(9);

  // test SetSqrtS()
  cout << "Testing SetSqrtS() ... ";
  myeB.SetSqrtS(200.0);
  gamma = myeB.GetGamma();
  Y0 = myeB.GetY0();
  if (fabs(gamma - 106.579) < tolerance && fabs(Y0 - 5.36201) < tolerance)
  {
    cout << "OK\n";
  }
  else
  {
    cout << "fail!"
         << " |gamma - 106.579| = " << fabs(gamma - 106.579)
         << " |Y0 - 5.36201| = " << fabs(Y0 - 5.36201) << endl;
  }

  // test SetGamma()
  cout << "Testing SetGamma() ... ";
  myeB.SetGamma(34.1052);
  sqrtS = myeB.GetSqrtS();
  Y0 = myeB.GetY0();
  if (fabs(sqrtS - 64.0) < tolerance && fabs(Y0 - 4.22238) < tolerance)
  {
    cout << "OK\n";
  }
  else
  {
    cout << "fail!\n";
  }

  // test SetY0
  cout << "Testing SetY0() ... ";
  myeB.SetY0(4.6688);
  sqrtS = myeB.GetSqrtS();
  gamma = myeB.GetGamma();
  if (fabs(sqrtS - 100.0) < tolerance && fabs(gamma - 53.2896) < tolerance)
  {
    cout << "OK\n";
  }
  else
  {
    cout << "fail!\n";
    cout << "    sqrtS - 100.0 = " << sqrtS - 100.0
         << "    Y0 - 4.6688 = " << gamma - 53.2896 << endl;
  }

  // 测试真空中磁场
  cout << "测试原点真空磁场: " << endl;
  double preCaleB[] = {68326.4, 215.59, 43.2741, 17.5784, 9.21409, 5.47008, 3.6497, 2.49731, 1.88101, 1.4337};

  myeB.SetSqrtS(200.0);
  for (int i = 0; i < 10; i++)
  {
    t = i * 5.0 / 10.0;
    myeB.SetSpaceTime(0.0, 0.0, 0.0, t);
    myeB.CalVaccumEB();
    cout << "t = " << myeB.GetT() << " eBy = " << myeB.eBy << " ... ";
    if (fabs(myeB.eBy - preCaleB[i]) / preCaleB[i] < 0.01)
    {
      cout << "OK";
    }
    else
    {
      cout << "fail";
    }
    cout << endl;
  }

  // 测试原点QGP响应
  cout << "测试原点QGP响应磁场" << endl;
  double preCalQGPeB[] = {68326.4, 2714.67, 1338.61, 871.988, 633.137, 485.838, 384.763, 310.534, 253.488, 208.269};

  for (int i = 0; i < 10; i++)
  {
    t = i * 5.0 / 10.0;
    myeB.SetSpaceTime(0.0, 0.0, 0.0, t);
    if (t <= myeB.GetTau0())
    {
      myeB.CalVaccumEB();
    }
    else
    {
      myeB.CalOriginQGPeB();
    }
    cout << "t = " << myeB.GetT() << " eBy = " << myeB.eBy << " ... ";
    if (fabs(myeB.eBy - preCalQGPeB[i]) / preCalQGPeB[i] < 0.01)
    {
      cout << "OK";
    }
    else
    {
      cout << "fail";
    }
    cout << endl;
  }

  myeB.SetSpaceTime_tau(0.25, 0.1);
  myeB.CalQGPeB();
  cout << "x = " << myeB.GetX() << " "
       << "y = " << myeB.GetY() << " "
       << "z = " << myeB.GetZ() << " "
       << "t = " << myeB.GetT() << endl
       << "eBy = " << myeB.eBy << endl;

  // 测试 SetTau0byCen()
  cout << "测试SetTau0byCen()函数" << endl;
  double tau0;
  myeB.SetNucleiType("Au");
  myeB.SetSqrtS(200.0);
  myeB.SetTau0byCen(35.0, 8.80);
  tau0 = myeB.GetTau0();
  cout << "Au-Au sqrts = 200 GeV centrality = 35% tau0 = " << tau0 << " ... ";
  if (fabs(tau0 - 0.162) < 0.001)
  {
    cout << "OK";
  }
  else
  {
    cout << "fail";
  }
  cout << endl;
  myeB.SetNucleiType("Pb");
  myeB.SetSqrtS(2760.0);
  myeB.SetTau0byCen(7.5, 4.31);
  tau0 = myeB.GetTau0();
  cout << "Pb-Pb sqrts = 2760 GeV centrality = 7.5% tau0 = " << tau0 << " ... ";
  if (fabs(tau0 - 0.095) < 0.001)
  {
    cout << "OK";
  }
  else
  {
    cout << "fail";
  }
  cout << endl;
  myeB.SetNucleiType("Cu");
  myeB.SetSqrtS(200.0);
  myeB.SetTau0byCen(15.0, 3.97);
  tau0 = myeB.GetTau0();
  cout << "Cu-Cu sqrts = 200 GeV centrality = 15% tau0 = " << tau0 << " ... ";
  if (fabs(tau0 - 0.165) < 0.001)
  {
    cout << "OK";
  }
  else
  {
    cout << "fail";
  }
  cout << endl;

  // 测试 cmefun
  cout << "测试cmefun()" << endl;
  myeB.SetMethod(0); // 0 for ellipsoid, 1 for disklike
  myeB.SetNucleiType("Au");
  myeB.SetSqrtS(200.0);
  myeB.SetTau0byCen(7.5, 4.05);
  myeB.SetNpm(577.6 / 2.0);
  myeB.SetLambda(0.2);
  myeB.CaleBy00();
  myeB.cmefun();
  if (fabs(myeB.app - 3.3e-8) < 0.1e-8 && fabs(myeB.apm - (-2.1e-9)) < 0.1e-9)
  {
    cout << "cmefun() 测试成功" << endl;
  }

  // test eB_Part_Int
  // int ndim = 4;
  // double xx[] = {0.4, 0.7, 0.2, 0.6};
  // int ncomp = 1;
  // double ff = 0.0;
  // myeB.flag = '+';
  // eB_Part_Int(&ndim, xx, &ncomp, &ff, (void *)&myeB);
  // std::cout << "eB_Part_Int = " << ff << std::endl;

  // test eB_Spec_Int
  // int ndim2 = 3;
  // double xx2[] = {0.2, 0.6, 0.4};
  // myeB.flag = '+';
  // eB_Spec_Int(&ndim2, xx2, &ncomp, &ff, (void *)&myeB);
  // std::cout << "eB_Spec_Int = " << ff << std::endl;

  return 0;
}
