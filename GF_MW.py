# Calculate MW from GF at 1-loop EW

import cmath as math

# Physical input parameters - couplings
aem = 1/137.036
GF = 1.1663787*10**(-5)

# Physical input parameters - masses
MZ = 91.1876
MZ2 = MZ**2
MT2 = 172.76**2
MB2 = 4.65**2
MH2 = 125.1**2

# shift of the electromagnetic coupling
DeltaAlpha = 0.05954


# other input
Pi = math.pi
Zeta2 = Pi**2/6
Zeta3 = 1.2020569


# 1-loop deltaR function as defined in eq (31)
def deltaR(MW2):
    deltar = ((1728*MW2**5 + 1664*MW2**4*MZ2 - 108*MT2**2*(2*MW2 - MZ2)*MZ2**2 -
               2*MW2**3*MZ2*(288*MB2 + 3329*MZ2) - 18*MZ2**3*(-6*MB2**2 + MH2**2 + MZ2**2) -
               9*MW2*MZ2**2*(24*MB2**2 - 4*MH2**2 - 6*MB2*MZ2 - 9*MH2*MZ2 + 27*MZ2**2) +
               MW2**2*MZ2*(-18*MH2**2 + 522*MB2*MZ2 - 81*MH2*MZ2 + 3527*MZ2**2) -
               6*MT2*(256*MW2**4 - 320*MW2**3*MZ2 + 73*MW2**2*MZ2**2 + 36*MB2*MZ2**3 -
                9*MW2*MZ2**2*(8*MB2 + MZ2)))/MW2**2 -
              (18*math.sqrt(-(MH2*(MH2 - 4*MW2)))*(MH2**2 - 4*MH2*MW2 + 12*MW2**2)*(2*MW2 - MZ2)*MZ2**2*
               math.atan(MH2/math.sqrt(-MH2**2 + 4*MH2*MW2)))/MW2**3 -
              (18*math.sqrt(-(MH2*(MH2 - 4*MW2)))*(MH2**2 - 4*MH2*MW2 + 12*MW2**2)*(2*MW2 - MZ2)*MZ2**2*
               math.atan((-MH2 + 2*MW2)/math.sqrt(-MH2**2 + 4*MH2*MW2)))/MW2**3 -
              (108*(MB2**2 + MT2**2 + MT2*MW2 - 2*MW2**2 + MB2*(-2*MT2 + MW2))*
               math.sqrt(-MB2**2 - (MT2 - MW2)**2 + 2*MB2*(MT2 + MW2))*(2*MW2 - MZ2)*MZ2**2*
               math.atan((-MB2 + MT2 - MW2)/math.sqrt(-MB2**2 - (MT2 - MW2)**2 + 2*MB2*(MT2 + MW2))))/MW2**3 +
              (108*(MB2**2 + MT2**2 + MT2*MW2 - 2*MW2**2 + MB2*(-2*MT2 + MW2))*
               math.sqrt(-MB2**2 - (MT2 - MW2)**2 + 2*MB2*(MT2 + MW2))*(2*MW2 - MZ2)*MZ2**2*
               math.atan((-MB2 + MT2 + MW2)/math.sqrt(-MB2**2 - (MT2 - MW2)**2 + 2*MB2*(MT2 + MW2))))/MW2**3 -
              72*math.sqrt(-1 + (4*MB2)/MZ2)*MZ2*(-8*MB2*MW2 + 11*MB2*MZ2 - 4*MW2*MZ2 + MZ2**2)*
              math.atan(1/math.sqrt(-1 + (4*MB2)/MZ2)) + 24*math.sqrt(-1 + (4*MT2)/MZ2)*
              (MT2*(64*MW2**2 - 80*MW2*MZ2 + 7*MZ2**2) + MZ2*(32*MW2**2 - 40*MW2*MZ2 + 17*MZ2**2))*
              math.atan(1/math.sqrt(-1 + (4*MT2)/MZ2)) + 36*math.sqrt(-1 + (4*MW2)/MZ2)*
              (-48*MW2**3 - 68*MW2**2*MZ2 + 16*MW2*MZ2**2 + MZ2**3)*math.atan(1/math.sqrt(-1 + (4*MW2)/MZ2)) +
              18*math.sqrt(-(MH2*(MH2 - 4*MZ2)))*(MH2**2 - 4*MH2*MZ2 + 12*MZ2**2)*
              math.atan(MH2/math.sqrt(-MH2**2 + 4*MH2*MZ2)) + 18*math.sqrt(-(MH2*(MH2 - 4*MZ2)))*
              (MH2**2 - 4*MH2*MZ2 + 12*MZ2**2)*math.atan((-MH2 + 2*MZ2)/math.sqrt(-MH2**2 + 4*MH2*MZ2)) +
              (18*MZ2*math.sqrt((4*MW2 - MZ2)*MZ2)*(-2*MW2 + MZ2)*(-48*MW2**3 - 68*MW2**2*MZ2 +
                16*MW2*MZ2**2 + MZ2**3)*math.atan(MZ2/math.sqrt(-MW2**2 - (-MW2 + MZ2)**2 + 2*MW2*(MW2 + MZ2))))/
              MW2**3 - (18*MZ2*math.sqrt((4*MW2 - MZ2)*MZ2)*(-2*MW2 + MZ2)*
                (-48*MW2**3 - 68*MW2**2*MZ2 + 16*MW2*MZ2**2 + MZ2**3)*
                math.atan((-2*MW2 + MZ2)/math.sqrt(-MW2**2 - (-MW2 + MZ2)**2 + 2*MW2*(MW2 + MZ2))))/MW2**3 +
              (54*(MB2 - MT2 + MW2)*(MB2**2 + MT2**2 + MT2*MW2 - 2*MW2**2 + MB2*(-2*MT2 + MW2))*
               (2*MW2 - MZ2)*MZ2**2*math.log(MB2/MT2))/MW2**3 -
              (9*MH2*(MH2**2 - 4*MH2*MW2 + 12*MW2**2)*(2*MW2 - MZ2)*MZ2**2*math.log(MH2/MW2))/MW2**3 +
              (36*MZ2**2*(3*MB2**2*(2*MT2 - MW2)*(2*MW2 - MZ2) + MT2*MW2**2*(-4*MW2 + MZ2) +
                MB2**3*(-6*MW2 + 3*MZ2) + MB2*(-3*MT2*MW2*(MW2 - 2*MZ2) + MW2**2*(4*MW2 - MZ2) +
                MT2**2*(-6*MW2 + 3*MZ2)))*math.log(MB2/MZ2))/((MB2 - MT2)*MW2**2) +
              9*MH2*(MH2**2 - 4*MH2*MZ2 + 12*MZ2**2 +
                     (2*(MW2 - MZ2)*MZ2*(MH2*MW2*(MW2 - 4*MZ2) + 12*MW2**2*MZ2 + MH2**2*(-MW2 + MZ2)))/
                     ((MH2 - MW2)*MW2**2))*math.log(MH2/MZ2) +
              (36*MZ2**2*(MB2**3*(6*MW2 - 3*MZ2) - 3*MB2**2*(2*MT2 - MW2)*(2*MW2 - MZ2) +
                         MT2*MW2**2*(4*MW2 - MZ2) + MB2*(MT2**2*(6*MW2 - 3*MZ2) + 3*MT2*MW2*(MW2 - 2*MZ2) +
                MW2**2*(-4*MW2 + MZ2)))*math.log(MT2/MZ2))/((MB2 - MT2)*MW2**2) +
              (9*MZ2**2*(2*MH2**3*MW2*(2*MW2 - MZ2) + 8*MH2**2*MW2**2*(-2*MW2 + MZ2) +
                MH2*(174*MW2**4 + 128*MW2**3*MZ2 - 130*MW2**2*MZ2**2 + 12*MW2*MZ2**3 + MZ2**4) -
                MW2*(144*MW2**4 + 152*MW2**3*MZ2 - 130*MW2**2*MZ2**2 + 12*MW2*MZ2**3 + MZ2**4))*
               math.log(MW2/MZ2))/(MW2**3*(-MH2 + MW2)))
    deltar = DeltaAlpha + deltar*aem/(Pi*(864*(MW2 - MZ2)**2*MZ2))
    return deltar

# solve for MW iteratively

# define start value for iteration
MW2it = 80**2

# accuracy for the iteration - start value and target value
acc = 1
acctarg = 0.00001

# iteratively solve for the MW fix point
while abs(acc) > acctarg:
    MW2old = MW2it
    MW2it = MZ2*(0.5 + math.sqrt(0.25 - aem*Pi/(math.sqrt(2)*GF*MZ2)*(1 + deltaR(MW2it))))
    acc = (1 - MW2it/MW2old).real


print(f"MW is {math.sqrt(MW2it).real}")