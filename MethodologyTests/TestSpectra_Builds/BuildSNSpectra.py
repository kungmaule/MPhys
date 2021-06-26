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
from DustAttenuation import *

##### METHODS ##################################################################



##### MAIN #####################################################################

if __name__=="__main__":

    DataSplineWindows = 'Parameter Files/Windows/RixWindows(Max2000).txt'

    ZArray = np.array([0.001, 0.002, 0.008, 0.014, 0.040])
    ZArraySolar = (ZArray/0.0142).round(2)

    # Iterate over all paths to generate the five models
    folderpath = r"C:\Users\hetms\Desktop\University\Year 5\MPhys\Code\sb99-v00-imf2.3-100myrs-models"
    filepaths = [os.path.join(folderpath, name) for name in os.listdir(folderpath)]
    ModelArray = np.empty(0)
    Stepper = 0
    print("Reading in Models...")
    for path in filepaths:
        Model = SpectrumModel.ReadInModelData(path, ZArraySolar[Stepper])
        ModelArray = np.append(ModelArray, Model)
        Stepper += 1

    # Read in example data for correct wavelength range
    SpectrumPath = f"Chi Square Fitting/Test_Spectra/test_spectrum0.dat"
    # Generate Model based upon test data
    StackModel = SpectrumModel.ReadInTestData(SpectrumPath)
    # Fit Spline to Model
    StackSpline = StackModel.DetermineContinuumSpline(DataSplineWindows)
    WavelengthGrid = StackSpline.Wavelengths

    # Define Metallicity Range
    ZRange = np.array([0.10, 0.25, 0.40])
    # Excludes start value as this already has a model
    #ZRange = ZRange[1:]
    SNArr = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15])

    for i in range(ZRange.shape[0]):
        Z = ZRange[i].round(2)
        if Z != 0.4:
            continue
        print(Z)
        for j in range(SNArr.shape[0]):
            SigToNoise = SNArr[j]
            if SigToNoise != 1:
                continue
            TestSpec = SpectrumModel.InterpolateModels(Z, ModelArray)
            # Add dust
            TestSpec.Errors = TestSpec.Fluxes / SigToNoise
            TestSpec.Fluxes = SpectrumModel.PerturbFluxes(TestSpec.Fluxes, TestSpec.Errors)
            TestSpec.ConvolveModel(FWHM=3.0, PixelRes=0.4)
            TestSpec.Fluxes = spectres(WavelengthGrid, TestSpec.Wavelengths, TestSpec.Fluxes)
            TestSpec.Errors = spectres(WavelengthGrid, TestSpec.Wavelengths, TestSpec.Errors)
            TestSpec.Wavelengths = WavelengthGrid

            Name = f'SNTest/TestSpectra/Z[{Z}]/test_spectrum{j}.dat'
            OutputTable = Table([TestSpec.Wavelengths, TestSpec.Fluxes, TestSpec.Errors], names=['wl', 'flam', 'flam_err'])
            OutputTable.write(Name, format='ascii', overwrite=True)
            print(f"Model of {Z}, SN={SigToNoise} generated.")
