#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "HIeB.hpp"
using namespace std;

int eBgen(const char Nuclei[], double sqrts, double b, double tau0, const char filename[], int type);

int main(int argc, char *argv[])
{
    double sqrts = 200.0;
    double b = 8.0;
    double tau0 = 0.154;
    int type = 0; //  0 for vaccum
    eBgen("Au", sqrts, b, tau0, "data/oriAu200b8Ai.dat", type);
    type = 1; // 1 for considering QGP
    eBgen("Au", sqrts, b, tau0, "data/oriAu200b8QGP.dat", type);
    
    sqrts = 2760.0;
    b = 8.0;
    tau0 = 0.108;
    type = 0;
    eBgen("Pb", sqrts, b, tau0, "data/oriPb2760b8Ai.dat", type);
    type = 1;
    eBgen("Pb", sqrts, b, tau0, "data/oriPb2760b8Ai.dat", type);

    /*
    cout << "Pb-Pb collisions, sqrts = 2760 GeV, magnetic field in vaccum generating ... ";
    myeB.SetNucleiType("Pb");
    myeB.SetSqrtS(2760.0);
    myeB.SetB(8.0);
    tau0 = 0.108;
    myeB.SetTau0(tau0);
    myeB.SetMethod(0); // 0 for ellipsoid, 1 for disklike

    ofstream output3("oriPb2760b8Ai.dat");
    t = 0.0;
    while (t < 5.0)
    {
        myeB.SetSpaceTime(0.0, 0.0, 0.0, t);
        myeB.CalVaccumEB();
        output3 << t << ", " << myeB.eBy << endl;
        t += 0.01;
    }
    output3.close();
    cout << "OK" << endl;

    cout << "Pb-Pb collisions, sqrts = 2760 GeV, magnetic field considering QGP generating ... ";
    myeB.SetNucleiType("Pb");
    myeB.SetSqrtS(2760.0);
    myeB.SetB(8.0);
    tau0 = 0.108;
    myeB.SetTau0(tau0);
    myeB.SetMethod(0); // 0 for ellipsoid, 1 for disklike

    ofstream output4("oriPb2760b8QGP.dat");
    t = 0.0;
    while (t < 5.0)
    {
        myeB.SetSpaceTime(0.0, 0.0, 0.0, t);
        if (t <= tau0)
        {
            myeB.CalVaccumEB();
            output4 << t << ", " << myeB.eBy << endl;
        }
        else
        {
            myeB.CalOriginQGPeB();
            output4 << t << ", " << myeB.eBy << endl;
        }
        t += 0.01;
    }
    output4.close();
    cout << "OK" << endl;
 */
    return 0;
}

int eBgen(const char Nuclei[], double sqrts, double b, double tau0, const char filename[], int type)
{
    HIeB myeB;

    myeB.SetNucleiType(Nuclei);
    myeB.SetSqrtS(sqrts);
    myeB.SetB(b);
    myeB.SetTau0(tau0);
    myeB.SetMethod(0); // 0 for ellipsoid, 1 for disklike

    cout << Nuclei << "-" << Nuclei << " collisions, sqrts = " << sqrts << "GeV,";
    if (type == 0) {
        cout << " magnetic field in vaccum generating ... ";
    } else {
        cout << " magnetic field considering QGP generating ... ";
    }
     

    ofstream output(filename);
    double t = 0.0;
    while (t < 5.0)
    {
        myeB.SetSpaceTime(0.0, 0.0, 0.0, t);
        if (type == 1 && t > tau0)
        { // 如果考虑QGP响应，并t>tau0时计算QGP中磁场
            myeB.CalOriginQGPeB();
        }
        else
        {
            myeB.CalVaccumEB();
        }
        output << t << ", " << myeB.eBy << endl;
        t += 0.01;
    }
    output.close();
    cout << "OK" << endl;

    return 0;
}