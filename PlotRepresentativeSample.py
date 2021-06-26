# Relevant Imports
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
from matplotlib import rc
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)

##### MAIN #####################################################################

if __name__=="__main__":

    # Read in FITS
    hdu = fits.open('StackingMethods/VANDWLS_ALL_DERIVED.fits')
    data = hdu[1].data
    hdr = hdu[0].header

    VlogMass = hdu[1].data['log10(M*)']
    VlogSFR = hdu[1].data['log10(SFR_UV_COR)']

    # Read in My Sample
    #Sample = Table.read('StackingMethods/params.sedparams', format='ascii')
    Sample = Table.read('StackingMethods/vandels_tstanton_dr3.dat', format='ascii')
    SlogMass = Sample['lmass']
    SlogSFR = Sample['lsfr']
    #SlogMass = Sample['lmass_fpp']
    #SlogSFR = Sample['lsfr_fpp']

    # Plot Figure
    plt.figure(figsize=(8,8))
    h = plt.hist2d(VlogMass, VlogSFR, bins=1000, label='Full VANDELS', cmap='gray_r')
    plt.scatter(SlogMass, SlogSFR, s=0.5, c='red', label='Sample')
    plt.xlim(7, 12)
    plt.ylim(-1, 3.5)
    plt.title('Plot of Sample over full VANDELS set')
    plt.xlabel('$\mathrm{\log(M_{*}/M_{\odot})}$')
    plt.ylabel('$\mathrm{\log(SFR/M_{\odot}yr^{-1})}$')
    plt.legend(loc='upper left')
    plt.colorbar(h[3], fraction=0.046, pad=0.04)
    plt.savefig('FINPLOT/RepresentativeSample.png')
    #plt.show()
