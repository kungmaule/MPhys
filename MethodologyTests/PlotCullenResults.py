# Relevant Imports
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib import rc, colors, cm
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)

##### MAIN #####################################################################

if __name__=="__main__":

    plt.figure(figsize=(8,8))
    #plt.xlim(0, 0.65)
    #plt.ylim(0, 0.65)
    #plt.plot([0, 0.65], [0, 0.65], lw=.5, color='black', alpha=.5, linestyle='dashed', label='1:1')

    # Read in Data
    Metallicity = np.zeros(7)
    LowZ = np.zeros(7)
    HighZ = np.zeros(7)
    lowerlimits = [False]*7
    upperlimits = [True] + [False]*6
    with open(f'Cullen2019/FCTestResults.txt', 'r') as file:
        for count, line in enumerate(file):
            Tokens = line.split(', ')
            LowZ[count] = np.log10(float(Tokens[2])) - np.log10(float(Tokens[1]))
            Metallicity[count] = np.log10(float(Tokens[2]))
            HighZ[count] = np.log10(float(Tokens[3])) - np.log10(float(Tokens[2]))

    ZError = np.vstack((LowZ, HighZ))
    F_Vals = np.array([-2.89, -2.78, -2.78, -2.70, -2.71, -2.56, -2.42])
    F_Errs = np.array([0.03, 0.03, 0.03, 0.03, 0.03, 0.06, 0.06])

    a, b = np.polyfit(F_Vals, Metallicity, 1)

    plt.xlim(-3.0, -2.2)
    plt.ylim(-3.0, -1.8)
    plt.errorbar(F_Vals, Metallicity, yerr=ZError, xerr=F_Errs, fmt='x', color='#D62839', capsize=2, linewidth=.5, label=f'C19', uplims=upperlimits, lolims=lowerlimits, xlolims=lowerlimits, xuplims=upperlimits)
    xlims = plt.gca().get_xlim()
    ylims = plt.gca().get_ylim()
    plt.plot(ylims, ylims, lw=.5, color='black', alpha=.5, linestyle='dashed', label='1:1')
    plt.plot(ylims, a*np.array(ylims) + b, color='#D62839', alpha=.5, label=f'C19 Fit: {a.round(2)}', lw=.5)
    plt.xlabel('C19 Metallicity / $\log (Z_{*})$')
    plt.ylabel('Recovered Metallicity / $\log (Z_{*})$')
    plt.legend(loc='upper left')
    plt.savefig('../FINPLOT/C19Comparison.png')
