#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "HIeB.hpp"
using namespace std;

#define LEN 9
int cme_ori(const char filename[], const char Nuclei[], double sqrts, double lambda, double cen[], double b[], double Npm2[]);

int main(int argc, char *argv[])
{
    //  Au-Au sqrts = 200 GeV
    double sqrts = 200.0;
    // double lambda = 0.2;
    double cen_l[LEN] = {0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0};
    double cen_u[LEN] = {5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0};
    double cen[LEN];
    for (int i = 0; i < LEN; i++)
    {
        cen[i] = (cen_l[i] + cen_u[i]) / 2, 0;
    }
    // b data from J. Phys. G: Nucl. Part. Phys. 35 (2008) 125106
    // double b_Au[LEN] = {2.30, 4.05, 5.75, 7.43, 8.80, 9.98, 11.04, 12.03, 12.96};
    // b data from PRC 79, 034909 (2009)
    double b_Au[LEN] = {2.21, 4.03, 5.70, 7.37, 8.73, 9.90, 11.0, 11.9, 12.8};
    double Npm2_Au[LEN] = {705.6, 577.6, 433.5, 291.5, 190.2, 118.5, 69.7, 37.8, 18.7};

    cme_ori("data/Au200GeV0.1.dat", "Au", sqrts, 0.1, cen, b_Au, Npm2_Au);
    cme_ori("data/Au200GeV0.2.dat", "Au", sqrts, 0.2, cen, b_Au, Npm2_Au);
    cme_ori("data/Au200GeV0.3.dat", "Au", sqrts, 0.3, cen, b_Au, Npm2_Au);

    // Pb-Pb sqrts = 2760 GeV
    sqrts = 2760.0;
    // b data inferd from PRC 88, 044909 (2013)
    double b_Pb[LEN] = {2.43, 4.31, 6.05, 7.81, 9.23, 10.47, 11.58, 12.58, 13.52};
    double Npm2_Pb[LEN] = {1601.0, 1294.0, 966.0, 649.0, 426.0, 261.0, 149.0, 76.0, 35.0};
    cme_ori("data/Pb2760GeV0.1.dat", "Pb", sqrts, 0.1, cen, b_Pb, Npm2_Pb);
    cme_ori("data/Pb2760GeV0.2.dat", "Pb", sqrts, 0.2, cen, b_Pb, Npm2_Pb);
    cme_ori("data/Pb2760GeV0.3.dat", "Pb", sqrts, 0.3, cen, b_Pb, Npm2_Pb);

    // Cu-Cu sqrts = 200 GeV
    sqrts = 200.0;
    // b data inferd from J. Phys. G: Nucl. Part. Phys. 35 (2008) 125106
    double b_Cu[LEN] = {1.75, 2.80, 3.97, 5.15, 6.10, 6.92, 7.68, 8.40, 9.10};
    double Npm2_Cu[LEN] = {184.8, 152.2, 117.4, 81.5, 55.4, 36.6, 23.3, 14.3, 8.4};
    cme_ori("data/Cu200GeV0.1.dat", "Cu", sqrts, 0.1, cen, b_Cu, Npm2_Cu);
    cme_ori("data/Cu200GeV0.2.dat", "Cu", sqrts, 0.2, cen, b_Cu, Npm2_Cu);
    cme_ori("data/Cu200GeV0.3.dat", "Cu", sqrts, 0.3, cen, b_Cu, Npm2_Cu);

    return 0;
}

int cme_ori(const char filename[], const char Nuclei[], double sqrts, double lambda, double cen[], double b[], double Npm2[])
{
    HIeB myeB;

    ofstream output(filename);
    myeB.SetMethod(1); // 0 for ellipsoid, 1 for disklike
    myeB.SetNucleiType(Nuclei);
    myeB.SetSqrtS(sqrts);
    myeB.SetLambda(lambda);
    output << "# app, apm, abs(apm)/app" << endl;
    for (int i = 0; i < LEN; i++)
    {
        myeB.SetTau0byCen(cen[i], b[i]);
        myeB.SetNpm(Npm2[i] / 2.0);
        myeB.CaleBy00();
        myeB.cmefun();
        output << myeB.app << ", " << myeB.apm << ", " << fabs(myeB.apm) / myeB.app << endl;
    }
    output.close();
    return 0;
}