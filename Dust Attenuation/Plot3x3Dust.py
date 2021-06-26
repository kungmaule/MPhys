import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.table import Table
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)
##### METHODS ##################################################################

def ReadInFile(Av, d, B):
    OutZ = np.empty(0)
    LowErr = np.empty(0)
    HighErr = np.empty(0)

    File = open(f'Data/ParameterTest/Av{Av}d{d}B{B}.txt', 'r')
    for line in File:
        Tokens = line.split(', ')
        OutZ = np.append(OutZ, float(Tokens[1]))
        LowErr = np.append(LowErr, float(Tokens[2]))
        HighErr = np.append(HighErr, float(Tokens[3]))

    Error = np.vstack((LowErr, HighErr))

    return OutZ, Error

def Residual(Value, TrueValue):
    return Value - TrueValue

##### MAIN #####################################################################

if __name__=="__main__":

    # Input Metallicities
    InZ = np.array([0.113, 0.156, 0.199, 0.242, 0.285, 0.328, 0.371, 0.414, 0.457, 0.5])

    rows, cols = 3, 3
    fig, ax = plt.subplots(rows, cols, sharex=True, sharey=True, figsize=(8,8))

    # Add plot data
    AvArr = [0, 2, 4]
    dArr = [-0.5, 0, 0.5]
    BArr = [1, 2, 4]
    c = ['#B8DE29', '#29AF7F', '#2D708E', '#482677']
    for row in range(rows):
        for col in range(cols):
            if row == 0:
                # Set out parameters
                d = 0
                B = 0
                Av = AvArr[col]
                # Read in Data
                OutZ, Error = ReadInFile(Av, d, B)
            elif row == 1:
                # Set out parameters
                d = dArr[col]
                B = 0
                Av = 1
                # Read in Data
                OutZ, Error = ReadInFile(Av, d, B)
            elif row == 2:
                # Set out parameters
                d = 0
                B = BArr[col]
                Av = 1
                # Read in Data
                OutZ, Error = ReadInFile(Av, d, B)

            # Plot
            label = '$A_{V}$= %.1f, $\delta$= %.1f, B= %.1f' % (Av, d, B)
            ax[row, col].errorbar(InZ, Residual(OutZ, InZ), yerr=Error, fmt='none', capsize=4, lw=.5, color=c[2])
            ax[row, col].scatter(InZ, Residual(OutZ, InZ), c=c[2], s=.4, label = label)
            ax[row, col].text(0.08, 0.4, f'Mean: %.4f, Variance: %.4f' % (np.mean(Residual(OutZ, InZ)), np.std(Residual(OutZ, InZ))),bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 5}, fontsize=9)
            ax[row, col].hlines(0.0, 0.06, 0.51, color='black', lw=.1)
            ax[row, col].legend(loc='lower left', prop={'size': 9})
            ax[row, col].set_ylim(-0.5, 0.5)
            ax[row, col].set_xlim(0.06, 0.51)

    ax[2,1].set_xlabel('Input Metallicity / $Z\\textsubscript{\(\odot\)}$')
    #ax2.set_ylabel('$A_{\lambda}$')
    fig.text(0.04, 0.5, 'Residual Metallicity / $Z\\textsubscript{\(\odot\)}$', va='center', rotation='vertical')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('../../FINPLOT/DustTest.png')
    plt.show()
