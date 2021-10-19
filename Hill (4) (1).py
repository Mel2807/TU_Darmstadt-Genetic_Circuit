import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Parameter import QscR_param, Membrane_Diff_param, Ahl_reg_param
from fractions import Fraction

# bas = 0.2, ktranskr = 2.04, ktransl = [2,5,4], Cn = 17, dmrna = 0.3301, dprot = 0.015, Kd_TF = 3.1, Kd_TFP = 0.53, nh = 1.7

def calc_epr_conc(nSig, basexpr, ktranskr, ktransl_rate, Cn, dmrna, dprot, Kd_TF, Kd_TFP, nh):
    '''
    calculates the concentration of the output-protein depending on the concentration of the signal molecule and the expression parameters
    @nSig: intracellular number of the signal molecule
    @basexpr: basal expression rate of the promoter for the output-protein
    @ktranskr: transkription rate of the output-protein
    @ktransl_rate: translation rate of the output-protein
    @Cn: copynumber of genetic circuit
    @dmrna: degradation rate of mRNA
    @dprot: degradation rate of output-protein
    @Kd_TF: Equilibrium constant of the dissociation of the transcription factor complex
    @Kd_TFP: Equilibrium constant of the dissociation of promotor and transcription factor complex
    @nh: hill coefficient
    '''
    # expr_max = (ktransl_rate * ktranskr * Cn) / (dmrna * dprot)
    # nSig = Fraction.from_float(nSig)
    # affinity_n = (Kd_TFP*Kd_TF*Cn*nSig.denominator**nh)
    # affinity_d = nSig.numerator**nh
    
    expr_max = (ktranskr * ktransl_rate * Cn) / (dmrna * dprot)


    affinity = (Kd_TFP * Kd_TF * Cn) / (nSig ** nh)
    c_out = expr_max * (basexpr + (1 - basexpr) * (nSig ** nh) / ((affinity) ** nh + nSig ** nh))
    return c_out

"""
complex = theta = (nsig)/(Kd_TF + nsig)

ODE
0 = K_TF_on * nsig^nh - K_TF_off * complex
complex = K_TF * nsig
"""


fig, ax = plt.subplots()
a = -2

nsig = np.logspace(14,17,100000)
nsig_cell = nsig * 6.5*10**(-16)


ktransl = np.linspace(2, 2, 1)
for k in range(len(ktransl)):
    nrun = np.copy(nsig_cell)

    # calculation of the protein output:
    for j in range(len(nrun)):
        nrun[j] = calc_epr_conc(nrun[j], QscR_param['basexpr'], QscR_param['ktranskr'], ktransl[k], QscR_param['Cn'], QscR_param['dmrna'], QscR_param['dprot'], QscR_param['Kd_TF'], QscR_param['Kd_TFP'], QscR_param['nh'])
    # nrun /= (6.5*10**(-10))*(6.2*10**23)*10**(-9)

    # save for plot
    ax.plot(nsig, nrun, label=ktransl[k], marker='.')


# PLOT THE OUTPUT:
# ax.set_ylim(0)
# ax.set_xlim(0)
ax.set_xlabel(r'#ahl_ext/L')
ax.set_ylabel(r'#protein/L')
ax.set_xscale('log')
ax.legend(loc="best", title="translation rates [1/min]", frameon=False)
# alternativ geht auch loc = 'best' für die Legendenposition
plt.show()