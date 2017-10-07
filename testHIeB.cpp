#include <iostream>
#include <iomanip>
#include <cmath>
#include "HIeB.hpp"

int main(int argc, char* argv[])
{
  double sqrtS, gamma, Y0, tolerance = 1e-3;
  HIeB myeB;

  std::setprecision(9);
  
  // test SetSqrtS()
  std::cout << "Testing SetSqrtS() ... ";
  myeB.SetSqrtS(200.0);
  gamma = myeB.GetGamma();
  Y0 = myeB.GetY0();
  if (fabs(gamma - 106.579) < tolerance && fabs(Y0 - 5.36201) < tolerance) {
    std::cout << "OK\n";
  } else {
    std::cout << "fail!" << " |gamma - 106.579| = " << fabs(gamma - 106.579)
	      << " |Y0 - 5.36201| = " << fabs(Y0 - 5.36201) << std::endl;
  }
  
  // test SetGamma()
  std::cout << "Testing SetGamma() ... ";
  myeB.SetGamma(34.1052);
  sqrtS = myeB.GetSqrtS();
  Y0 = myeB.GetY0();
  if (fabs(sqrtS - 64.0) < tolerance && fabs(Y0 - 4.22238) < tolerance) {
    std::cout << "OK\n";
  } else {
    std::cout << "fail!\n";
  }

  // test SetY0
  std::cout << "Testing SetY0() ... ";
  myeB.SetY0(4.6688);
  sqrtS = myeB.GetSqrtS();
  gamma = myeB.GetGamma();
  if (fabs(sqrtS - 100.0) < tolerance && fabs(gamma - 53.2896) < tolerance) {
    std::cout << "OK\n";
  } else {
    std::cout << "fail!\n";
    std::cout << "    sqrtS - 100.0 = " << sqrtS - 100.0
	      << "    Y0 - 4.6688 = " << gamma - 53.2896 << std::endl;
  }

  /*
  int i;
  double sign;
  myeB.x = 0.0;
  myeB.y = 0.0;
  myeB.z = 0.0;
  myeB.xp = -6.0;
  myeB.yp = 0.0;
  myeB.SetSqrtS(200.0);
  sign = 1.0;

  std::cout << "gamma = " << myeB.GetGamma() << "\n";
  std::cout << "sigma = " << myeB.sigma << "\n";
  std::cout << "sigmaChi = " << myeB.sigmaChi << "\n";
  std::cout << "v = " << myeB.GetV() << "\n";
  for (i = 0; i < 100; i++) {
    myeB.t = i/100.0 * 1.0;
    myeB.pointEB_Wangqun(sign);
    std::cout << "t = " << myeB.t << " eBy = " << myeB.GetPointEB(2) << "\n";
  }
  */
}
