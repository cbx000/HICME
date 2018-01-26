# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from numpy import genfromtxt, savetxt
from scipy.optimize import minimize_scalar

# read exp gamma data
Au62opp = genfromtxt('../data/exp/Au62GeVopp.csv', delimiter=',')
Au62same = genfromtxt('../data/exp/Au62GeVsame.csv', delimiter=',')
Au200opp = genfromtxt('../data/exp/STARopp.csv', delimiter=',')
Au200same = genfromtxt('../data/exp/STARsame.csv', delimiter=',')
Pb2760opp = genfromtxt('../data/exp/ALICEopp.csv', delimiter=',')
Pb2760same = genfromtxt('../data/exp/ALICEsame.csv', delimiter=',')
Cu62opp = genfromtxt('../data/exp/Cu62GeVopp.csv', delimiter=',')
Cu62same = genfromtxt('../data/exp/Cu62GeVsame.csv', delimiter=',')
Cu200opp = genfromtxt('../data/exp/Cu200GeVopp.csv', delimiter=',')
Cu200same = genfromtxt('../data/exp/Cu200GeVsame.csv', delimiter=',')

# read exp delta data
Au200deltasame = genfromtxt('../data/exp/Au200GeV_delta_same.txt', delimiter=',')
Au200deltaopp = genfromtxt('../data/exp/Au200GeV_delta_opp.txt', delimiter=',')
Pb2760deltasame = genfromtxt('../data/exp/Pb2760GeV_delta_same.txt', delimiter=',')
Pb2760deltaopp = genfromtxt('../data/exp/Pb2760GeV_delta_opp.txt', delimiter=',')

# read exp v2 data
Au200v2 = genfromtxt('../data/v2/RHICv2.txt', delimiter=',', comments='#')
Pb2760v2 = genfromtxt('../data/v2/LHCv2.txt', delimiter=',', comments='#')

# evaluate H
kappa = 1.5
Au200Hsame = (kappa * Au200v2[:,1] * Au200deltasame[1:,1] - Au200same[1:,1])\
    /(1 + kappa*Au200v2[:,1])
Au200Hopp = (kappa * Au200v2[:,1] * Au200deltaopp[1:,1] - Au200opp[1:,1])\
    /(1 + kappa*Au200v2[:,1])
Pb2760Hsame = (kappa * Pb2760v2[:,1] * Pb2760deltasame[:,1] - Pb2760same[:,1])\
    /(1 + kappa*Pb2760v2[:,1])
Pb2760Hopp = (kappa * Pb2760v2[:,1] * Pb2760deltaopp[:,1] - Pb2760opp[:,1])\
    /(1 + kappa*Pb2760v2[:,1])
Au200Hdiff = Au200Hsame - Au200Hopp
Pb2760Hdiff = Pb2760Hsame - Pb2760Hopp

savetxt('Au200Hsame.txt', Au200Hsame)
savetxt('Au200Hopp.txt', Au200Hopp)
savetxt('Pb2760Hsame.txt', Pb2760Hsame)
savetxt('Pb2760Hopp.txt', Pb2760Hopp)

# read theory data
Au20001 = genfromtxt('../data/Au200GeV0.1.dat', delimiter=',')
Au20002 = genfromtxt('../data/Au200GeV0.2.dat', delimiter=',')
Au20003 = genfromtxt('../data/Au200GeV0.3.dat', delimiter=',')
Pb276001 = genfromtxt('../data/Pb2760GeV0.1.dat', delimiter=',')
Pb276002 = genfromtxt('../data/Pb2760GeV0.2.dat', delimiter=',')
Pb276003 = genfromtxt('../data/Pb2760GeV0.3.dat', delimiter=',')

Au20001Diff = Au20001[1:8,0] - Au20001[1:8,1]
Au20002Diff = Au20002[1:8,0] - Au20002[1:8,1]
Au20003Diff = Au20003[1:8,0] - Au20003[1:8,1]
Pb276001Diff = Pb276001[:8,0] - Pb276001[:8,1]
Pb276002Diff = Pb276002[:8,0] - Pb276002[:8,1]
Pb276003Diff = Pb276003[:8,0] - Pb276003[:8,1]


def adjustfun(Hdiff, theDiff):
    def objfun(alpha):
        temp = (Hdiff - theDiff*alpha)**2/Hdiff
        return sum(temp)
    res = minimize_scalar(objfun, method='brent')
    print("alpha =", res.x, "chi2 =", res.fun)
    return res


# H diff result
print("Au 200GeV lambda/R = 0.1, 0.2, 0.3")
adjustfun(Au200Hdiff[:6], Au20001Diff[:6])
adjustfun(Au200Hdiff[:6], Au20002Diff[:6])
adjustfun(Au200Hdiff[:6], Au20003Diff[:6])

print("Pb 2760GeV lambda/R = 0.1, 0.2, 0.3")
adjustfun(Pb2760Hdiff[:7], Pb276001Diff[:7])
adjustfun(Pb2760Hdiff[:7], Pb276002Diff[:7])
adjustfun(Pb2760Hdiff[:7], Pb276003Diff[:7])
