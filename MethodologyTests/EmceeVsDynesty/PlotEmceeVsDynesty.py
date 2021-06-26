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

    """
    # Read in Emcee Data
    EmceeData = np.zeros(shape=(3,20))
    with open('EmceeResults.txt', 'r') as file:
        for count, line in enumerate(file):
            Tokens = line.split(', ')
            print(Tokens)
            EmceeData[1][count] = np.log10(float(Tokens[1]))
            EmceeData[0][count] = np.log10(float(Tokens[1])) - np.log10(float(Tokens[1])-float(Tokens[2]))
            EmceeData[2][count] = np.log10(float(Tokens[3])-float(Tokens[1])) - np.log10(float(Tokens[1]))
        EmceeError = np.vstack((EmceeData[0],EmceeData[2]))
        EmceeValues = EmceeData[1]

    # Read in Dynesty Data
    DynestyData = np.zeros(shape=(3,20))
    with open('DynestyResults.txt', 'r') as file:
        for count, line in enumerate(file):
            Tokens = line.split(', ')
            DynestyData[1][count] = (float(Tokens[2]))
            DynestyData[0][count] = (float(Tokens[2])) - (float(Tokens[1]))
            DynestyData[2][count] = (float(Tokens[3])) - (float(Tokens[2]))
        DynestyError = np.vstack((DynestyData[0], DynestyData[2]))
        DynestyValues = DynestyData[1]

    """
    # Read in Emcee Data
    EmceeData = np.zeros(shape=(3,20))
    with open('EmceeResults.txt', 'r') as file:
        for count, line in enumerate(file):
            Tokens = line.split(', ')
            EmceeData[1][count] = float(Tokens[1]) / 0.0142
            EmceeData[0][count] = float(Tokens[2]) / 0.0142
            EmceeData[2][count] = float(Tokens[3]) / 0.0142
        EmceeError = np.vstack((EmceeData[0],EmceeData[2]))
        EmceeValues = EmceeData[1]

    # Read in Dynesty Data
    DynestyData = np.zeros(shape=(3,20))
    with open('DynestyResults.txt', 'r') as file:
        for count, line in enumerate(file):
            Tokens = line.split(', ')
            DynestyData[1][count] = (float(Tokens[2]) / 0.0142)
            DynestyData[0][count] = (float(Tokens[2]) - float(Tokens[1])) / 0.0142
            DynestyData[2][count] = (float(Tokens[3]) - float(Tokens[2])) / 0.0142
        DynestyError = np.vstack((DynestyData[0], DynestyData[2]))
        DynestyValues = DynestyData[1]

    TrueValues = np.linspace(0.07, 0.5, 21)
    TrueValues = TrueValues[1:] #* 0.0142
    #TrueValues = np.log10(TrueValues)

    aDynesty, bDynesty = np.polyfit(TrueValues, DynestyValues, 1)
    aEmcee, bEmcee = np.polyfit(TrueValues, EmceeValues, 1)

    # Plot
    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True, figsize=(6,8), gridspec_kw={'height_ratios': [2,1]})
    plt.xlim(0, 0.65)
    ax1.set_ylim(0, 0.65)
    ax2.set_ylim(-.25, .25)
    ax1.plot([0, 0.65], [0, 0.65], lw=.5, color='black', alpha=.5, linestyle='dashed', label='1:1')
    #ax1.plot(TrueValues, TrueValues, lw=.5, color='black', alpha=.5, linestyle='dashed', label='1:1')
    # Plot Emcee Results
    ax1.errorbar(TrueValues, EmceeValues, yerr=EmceeError, fmt='x', color='#D62839', capsize=2, linewidth=.5, label=f'Emcee')
    ax1.plot(TrueValues, aEmcee*TrueValues + bEmcee, color='#D62839', alpha=.5, label=f'Emcee Fit: {aEmcee.round(2)}', lw=.5)

    # Plot Dynesty Results
    ax1.errorbar(TrueValues, DynestyValues, yerr=DynestyError, fmt='x', color='#00916E', capsize=2, linewidth=.5, label=f'Dynesty')
    ax1.plot(TrueValues, aDynesty*TrueValues + bDynesty, color='#00916E', alpha=.5, label=f'Dynesty Fit: {aDynesty.round(2)}', lw=.5)

    # Plot Residue
    ax2.errorbar(TrueValues, EmceeValues-TrueValues, yerr=EmceeError, fmt='x', color='#D62839', capsize=2, linewidth=.5, label=f'Emcee')
    ax2.errorbar(TrueValues, DynestyValues-TrueValues, yerr=DynestyError, fmt='x', color='#00916E', capsize=2, linewidth=.5, label=f'Dynesty')
    ax2.axhline(0.0, lw=.25, color='black')

    ax2.set_xlabel('True Metallicity / $Z_{\odot}$')
    ax1.set_ylabel('Recovered Metallicity / $ Z_{\odot}$')
    ax2.set_ylabel('Residual from 1:1')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper left')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('../../FINPLOT/EmceeVsDynesty.png')
