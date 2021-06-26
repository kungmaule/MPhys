import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from astropy.table import Table
from matplotlib import colors, cm
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)
##### METHODS ##################################################################

def ReadInFile(Z):
    OutZ = np.empty(0)
    LowErr = np.empty(0)
    HighErr = np.empty(0)

    File = open(f'Data/Z{Z}.txt', 'r')
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

    CArr = np.empty(0, dtype='str')
    cmap = cm.get_cmap('plasma', 9)
    for i in range(1, cmap.N-1):
        rgba = cmap(i)
        col = colors.rgb2hex(rgba)
        CArr = np.append(CArr, col)

    # Input Metallicities
    InZ = np.array([0.1, 0.25, 0.4])
    SNArr = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15])
    #SNArr = np.log10(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15]))

    rows, cols = 2, 1
    fig, ax = plt.subplots(rows, cols, sharex=True, sharey=False, figsize=(8,8))
    c = [CArr[0], CArr[3], CArr[6]]
    # Generate Top Plot
    for i in range(InZ.shape[0]):
        Z = InZ[i].round(2)
        col = c[i]
        OutZ, Error = ReadInFile(Z)
        Res = np.array([abs((ZDot-Z)/Z) for ZDot in OutZ])
        #SNArr = np.log10(SNArr)

        # Top Plot
        ax[0].plot(SNArr, OutZ, lw=.5, color=col)
        ax[0].errorbar(SNArr, OutZ, yerr=Error, fmt='none', capsize=4, lw=.5, color=col)
        ax[0].scatter(SNArr, OutZ , c=col, s=.6, label='$Z$/$Z\\textsubscript{\(\odot\)}$ = %.2f' % Z)
        ax[0].legend(loc='upper right')

        # Lower Plot
        def GuessFit(x, a, b, c):
            return a + b/x
            #return np.exp(-b*x) + c
            #return a*np.power(x,2) + b*x + c

        SNFit = np.linspace(SNArr[0], SNArr[-1], 100)
        popt, pcov = curve_fit(GuessFit, SNArr, Res)
        ResFit = GuessFit(SNFit, *popt)
        print(ResFit, SNFit)
        print('\n', SNFit[32], ResFit[32])
        ax[1].plot(SNFit, ResFit, lw=.5, color=col)
        #ax[1].plot(SNArr, Res, color=col, lw=.5)
        ax[1].scatter(SNArr, Res , c=col, s=.6, label='$Z$/$Z\\textsubscript{\(\odot\)}$ = %.2f' % Z)
        ax[1].legend(loc='upper right')
        ax[1].set_ylim(0.0, 6.5)
        ax[1].set_yticks([0, 1, 2, 3, 4, 5, 6])

    ax[1].set_xlabel('S/N Ratio')
    ax[0].set_ylabel('$Z$ / $Z\\textsubscript{\(\odot\)}$')
    ax[1].set_ylabel('$|Z-Z_{True}|$ / $Z_{True}$')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlim(0, 15)
    #plt.xlim(0, np.log10(15))
    plt.savefig('../../FINPLOT/SNPlot.png')
    plt.show()
