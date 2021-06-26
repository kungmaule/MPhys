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
    plt.xlim(0, 0.65)
    plt.ylim(0, 0.65)
    plt.plot([0, 0.65], [0, 0.65], lw=.5, color='black', alpha=.5, linestyle='dashed', label='1:1')

    Name = 'S99'
    # Read in Data
    Metallicity = np.zeros(20)
    LowZ = np.zeros(20)
    HighZ = np.zeros(20)
    with open(f'{Name}.txt', 'r') as file:
        for count, line in enumerate(file):
            Tokens = line.split(', ')
            LowZ[count] = float(Tokens[2]) - float(Tokens[1])
            Metallicity[count] = float(Tokens[2])
            HighZ[count] = float(Tokens[3]) - float(Tokens[2])

    Error = np.vstack((LowZ, HighZ))
    ZRange = np.linspace(0.07, 0.5, 21)
    # Excludes start value as this already has a model
    InZ = ZRange[1:]

    a, b = np.polyfit(InZ, Metallicity, 1)

    plt.errorbar(InZ, Metallicity, yerr=Error, fmt='x', color='#D62839', capsize=2, linewidth=.5, label=f'S99')
    plt.plot(InZ, a*InZ + b, color='#D62839', alpha=.5, label=f'S99 Fit: {a.round(2)}', lw=.5)

    Name = 'BPASS'
    # Read in Data
    Metallicity = np.zeros(20)
    LowZ = np.zeros(20)
    HighZ = np.zeros(20)
    with open(f'{Name}.txt', 'r') as file:
        for count, line in enumerate(file):
            Tokens = line.split(', ')
            LowZ[count] = float(Tokens[2]) - float(Tokens[1])
            Metallicity[count] = float(Tokens[2])
            HighZ[count] = float(Tokens[3]) - float(Tokens[2])

    a, b = np.polyfit(InZ, Metallicity, 1)

    plt.errorbar(InZ, Metallicity, yerr=Error, fmt='x', color='#00916E', capsize=2, linewidth=.5, label=f'BPASS: Grad {a.round(2)}')
    plt.plot(InZ, a*InZ + b, color='#00916E', alpha=.5, label=f'BPASS Fit: {a.round(2)}', lw=.5)
    plt.xlabel('True Metallicity / $Z_{\odot}$')
    plt.ylabel('Recovered Metallicity / $Z_{\odot}$')
    plt.legend(loc='upper left')
    plt.savefig('../../FINPLOT/S99vsBPASS.png')
