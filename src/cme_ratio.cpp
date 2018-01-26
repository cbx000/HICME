#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "HIeB.hpp"
using namespace std;

int cme_ratio(const char filename[], const char Nuclei[], double sqrts, double lambda);

int main(int argc, char *argv[])
{

    cme_ratio("data/ratioAu200GeV0.1.dat", "Au", 200.0, 0.1);
    cme_ratio("data/ratioAu200GeV0.2.dat", "Au", 200.0, 0.2);
    cme_ratio("data/ratioAu200GeV0.3.dat", "Au", 200.0, 0.3);
    return 0;
}

int cme_ratio(const char filename[], const char Nuclei[], double sqrts, double lambda)
{
    HIeB myeB;
    ofstream output(filename);
    myeB.SetMethod(1);  // 0 for ellipsoid, 1 for disklike
    myeB.SetNucleiType(Nuclei);
    myeB.SetSqrtS(sqrts);
    myeB.SetLambda(lambda);
    double bmin = 0.0;
    double R = myeB.GetR();
    double bmax = 2.0 * R;
    int N = 500;
    double b;
    output << "# b/R, abs(apm)/app" << endl;
    for (int i = 0; i < N; i++)
    {
        b = bmin + (bmax - bmin)*i/N;
        myeB.SetB(b);
        myeB.SetTau0(0.1); // tau0的值也无所谓，对磁场的积分会抵消掉
        myeB.SetNpm(2.0);// Npm的值无所谓，因为这里计算的是比值，最终会消去
        myeB.CaleBy00();
        myeB.cmefun();
        output << b/R << ", " << fabs(myeB.apm) / myeB.app << endl;
    }
    return 0;
}
