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
  myeB.SetSpaceTime(0.0, 0.0, 0.0, 0.1);
  myeB.SetSqrtS(200.0);
  
  myeB.CalVaccumEB();
  std::cout << "x = " << myeB.GetX() << std::endl
  	    << "y = " << myeB.GetY() << std::endl
  	    << "z = " << myeB.GetZ() << std::endl
  	    << "t = " << myeB.GetT() << std::endl
  	    << "eBy = " << myeB.eBy << std::endl;
  

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
