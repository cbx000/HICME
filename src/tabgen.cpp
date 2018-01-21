#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "HIeB.hpp"
using namespace std;
#define LEN 9

int procfun(char *filename, char *Nuclei, double sqrts, char *note,
            double cen_l[], double cen_u[], double cen[], double b[], int len);

int main(int argc, char *argv[])
{
    // Table 1. Au-Au sqrts = 200 GeV
    char filename[] = "data/tab1:Au200GeV.dat";
    char Nuclei[] = "Au";
    double sqrts = 200.0;
    char note[] = "%Table 1 Centrality dependence of Qs2, tau0, eBy0\n"
                  "% for Au-Au collisions, sqrts = 200 GeV\n"
                  "% Centrality, b, Qs2, tau0, eBy0";
    double cen_l[LEN] = {0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0};
    double cen_u[LEN] = {5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0};
    double cen[LEN];
    for (int i = 0; i < LEN; i++)
    {
        cen[i] = (cen_l[i] + cen_u[i]) / 2, 0;
    }
    // b data from J. Phys. G: Nucl. Part. Phys. 35 (2008) 125106
    double b_Au[LEN] = {2.30, 4.05, 5.75, 7.43, 8.80, 9.98, 11.04, 12.03, 12.96};
    procfun(filename, Nuclei, sqrts, note, cen_l, cen_u, cen, b_Au, LEN - 2);

    // Table 2. Pb-Pb sqrts = 2760 GeV
    char filename2[] = "data/tab2:Pb2760GeV.dat";
    char Nuclei2[] = "Pb";
    sqrts = 2760.0;
    char note2[] = "%Table 2 Centrality dependence of Qs2, tau0, eBy0\n"
                   "% for Pb-Pb collisions, sqrts = 2760 GeV\n"
                   "% Centrality, b, Qs2, tau0, eBy0";
    // b data inferd from PRC 88, 044909 (2013)
    double b_Pb[LEN] = {2.43, 4.31, 6.05, 7.81, 9.23, 10.47, 11.58, 12.58, 13.52};
    procfun(filename2, Nuclei2, sqrts, note2, cen_l, cen_u, cen, b_Pb, LEN - 2);

    // Table 3. Cu-Cu sqrts = 200 GeV
    char filename3[] = "data/tab3:Cu200GeV.dat";
    char Nuclei3[] = "Cu";
    sqrts = 200.0;
    char note3[] = "%Table 3 Centrality dependence of Qs2, tau0, eBy0\n"
                   "% for Cu-Cu collisions, sqrts = 200 GeV\n"
                   "% Centrality, b, Qs2, tau0, eBy0";
    // b data inferd from J. Phys. G: Nucl. Part. Phys. 35 (2008) 125106
    double b_Cu[LEN] = {1.75, 2.80, 3.97, 5.15, 6.10, 6.92, 7.68, 8.40, 9.10};
    procfun(filename3, Nuclei3, sqrts, note3, cen_l, cen_u, cen, b_Cu, LEN - 2);

    return 0;
}

int procfun(char *filename, char *Nuclei, double sqrts, char *note,
            double cen_l[], double cen_u[], double cen[], double b[], int len)
{
    HIeB myeB;
    double eBy0 = 0.0;
    double Qs2;

    ofstream out(filename);
    myeB.SetNucleiType(Nuclei);
    myeB.SetSqrtS(sqrts);
    out << note << endl;

    myeB.SetMethod(0); // ellipsoid

    for (int i = 0; i < len; i++)
    {
        myeB.SetB(b[i]);
        Qs2 = myeB.SetTau0byCen(cen[i], b[i]);
        myeB.CaleBy00();
        eBy0 = myeB.GeteBy00();
        out << (int)cen_l[i] << "--" << (int)cen_u[i] << "\\% & " << fixed << setprecision(2) << myeB.GetB() << " & " << fixed << setprecision(2) << Qs2 << " & " << setprecision(3) << myeB.GetTau0() << " & " << fixed << setprecision(1) << eBy0 << " \\\\" << endl;
    }
    out.close();
    return 0;
}