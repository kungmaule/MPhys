# Imports from libraries
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)
import os
import sys
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
from spectres import spectres
# Imports from Classes
from SpectrumModel import SpectrumModel, Spline, GenerateColours

##### METHODS ##################################################################

def A_lambda_mod(wl, Av, delta, B):
    RvCalz = 4.05 # Cullen 2019
    RvMod = RvCalz/((RvCalz+1)*((4400/5500)**delta)-RvCalz) # Salim 2018
    k = k_Calz(wl, RvCalz)
    k_mod = (k * (RvMod/RvCalz) * np.power(wl/5500, delta)) + D_lambda(wl, B)
    return (Av/RvMod)*k_mod

def D_lambda(wl, B):
    dwl = 305/10000
    wlc = 2175/10000
    wl = wl/10000
    d_lam = (B*np.power(wl*dwl,2)) / ((np.power(wl**2-wlc**2,2)) + np.power(wlc*dwl,2))
    return d_lam

def k_Calz(wl, RvCalz):
    if 1200 <= wl <= 6300:
        wl = wl * 1e-4
        return 2.659*(-2.156 + (1.509/wl) - (0.19/(wl**2)) + (0.011/(wl**3))) + RvCalz
    # Pulled from Calzetti's website ? May be wrong?
    # Can extrapolate back but shouldn't need to


##### MAIN #####################################################################

if __name__=="__main__":

    Wavelengths = np.arange(1200, 6300, 1)
    c = GenerateColours(Seq=True)

    AvRange = [0, 1, 2, 4]
    SlRange = [-0.5, 0.0, 0.5]
    BRange = [0, 1, 2, 4]
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(10,10))
    for i in range(len(AvRange)):
        Av = AvRange[i]
        B = 0.
        delta = 0.0
        A_lambda = [A_lambda_mod(wl, Av, delta, B) for wl in Wavelengths]
        ax1.plot(Wavelengths, A_lambda, color=c[i], lw=.5, label= '$A_{V}$= %.1f, $\delta$= %.1f, B= %.1f' % (Av, delta, B))
    for i in range(len(SlRange)):
        Av = 1.0
        B = 0.
        delta = SlRange[i]
        A_lambda = [A_lambda_mod(wl, Av, delta, B) for wl in Wavelengths]
        ax2.plot(Wavelengths, A_lambda, color=c[i], lw=.5, label= '$A_{V}$= %.1f, $\delta$= %.1f, B= %.1f' % (Av, delta, B))
    for i in range(len(BRange)):
        Av = 1.0
        B = BRange[i]
        delta = 0.0
        A_lambda = [A_lambda_mod(wl, Av, delta, B) for wl in Wavelengths]
        ax3.plot(Wavelengths, A_lambda, color=c[i], lw=.5, label= '$A_{V}$= %.1f, $\delta$= %.1f, B= %.1f' % (Av, delta, B))

    ax1.legend(loc='upper right')
    ax2.legend(loc='upper right')
    ax3.legend(loc='upper right')
    plt.xlabel('Wavelength / Å')
    #ax2.set_ylabel('$A_{\lambda}$')
    fig.text(0.04, 0.5, '$A_{\lambda}$', va='center', rotation='vertical')
    plt.xlim(1200, 5000)
    plt.subplots_adjust(wspace=0, hspace=0)
    #plt.ylim(0, 8)
    plt.savefig('../FINPLOT/DustAttenuationGraph.png')
    plt.show()
