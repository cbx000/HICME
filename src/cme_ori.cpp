#include <iostream>
#include <iomanip>
#include <cmath>
#include "HIeB.hpp"
using namespace std;

#define LEN 9

int main(int argc, char *argv[])
{
    HIeB myeB;
    //  Au-Au sqrts = 200 GeV
    double cen_l[LEN] = {0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0};
    double cen_u[LEN] = {5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0};
    double cen[LEN];
    for (int i = 0; i < LEN; i++)
    {
        cen[i] = (cen_l[i] + cen_u[i]) / 2, 0;
    }
    // b data from J. Phys. G: Nucl. Part. Phys. 35 (2008) 125106
    double b_Au[LEN] = {2.30, 4.05, 5.75, 7.43, 8.80, 9.98, 11.04, 12.03, 12.96};
    double Npm2_Au[LEN] = {705.6, 577.6, 433.5, 291.5, 190.2, 118.5, 69.7, 37.8, 18.7};

    myeB.SetMethod(0); // 0 for ellipsoid, 1 for disklike
    myeB.SetNucleiType("Au");
    myeB.SetSqrtS(200.0);
    myeB.SetLambda(0.2);
    printf("# app, apm, abs(apm)/app\n");
    for (int i = 0; i < LEN; i++)
    {
        myeB.SetTau0byCen(cen[i], b_Au[i]);
        myeB.SetNpm(Npm2_Au[i] / 2.0);
        myeB.CaleBy00();
        myeB.cmefun();
        printf("%g, %g, %g\n", myeB.app, myeB.apm, fabs(myeB.apm)/myeB.app);
    }

    return 0;
}