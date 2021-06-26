# Imports from libraries
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
import os
import sys
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
from spectres import spectres
# Imports from Classes
from SpectrumModel import SpectrumModel, Spline


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

##### MAIN #####################################################################

if __name__=="__main__":

    # Iterate over all paths to read in the five models
    folderpath = r"C:\Users\hetms\Desktop\University\Year 5\MPhys\Code\Dust Attenuation\TestSpectra\NonDusty"
    filepaths = [os.path.join(folderpath, name) for name in os.listdir(folderpath)]
    ModelArray = np.empty(0)
    Stepper = 0
    print("Reading in Data...")
    for path in filepaths:
        ModelArray = np.append(ModelArray, SpectrumModel.ReadInTestData(path))

    # Define Wavelength Range
    WavelengthGrid = ModelArray[0].Wavelengths

    # Generate Dust Attenuation Curve
    Av = 4
    delta = 0
    B = 0
    A_lambda = np.array([A_lambda_mod(wl, Av, delta, B) for wl in WavelengthGrid])

    # Add dust to the indiviudual modles

    for i in range(ModelArray.shape[0]):
        Model = ModelArray[i]
        Model.Fluxes = Model.Fluxes * np.exp(-0.4*A_lambda)
        Model.Errors = Model.Errors * np.exp(-0.4*A_lambda)

        Name = f'Dust Attenuation/TestSpectra/Av{int(Av)}d{delta}B{int(B)}/test_spectrum{i}.dat'
        OutputTable = Table([Model.Wavelengths, Model.Fluxes, Model.Errors], names=['wl', 'flam', 'flam_err'])
        OutputTable.write(Name, format='ascii', overwrite=True)
        print(f"Dust added to Model {i}.")
