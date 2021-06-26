# Relevant Imports
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib import rc
import dynesty as dyn
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)
from SpectrumModel import SpectrumModel, Spline, GenerateColours
##### MAIN #####################################################################

if __name__=="__main__":

    # Define Necessary Constants
    ModelSplineWindows = 'Parameter Files/Windows/RixWindows(Max2000).txt'
    DataSplineWindows = 'Parameter Files/Windows/RixWindows(Max2000).txt'
    ChiSquareWindows = 'Parameter Files/Windows/SteidelWindows.txt'

    # Generate Model Array
    ZArray = np.array([0.001, 0.002, 0.008, 0.014, 0.040])
    FWHM, Res = 3.0, 0.4 # 0.4Å Resolution
    folderpath = r"./S99"
    ResultPath = 'S99'

    #ZArraySolar = (ZArray/0.0142).round(2)
    filepaths = [os.path.join(folderpath, name) for name in os.listdir(folderpath)]
    ModelArray = np.empty(0)
    ModelWavelengths = np.arange(1000, 2000, 1)
    for count, path in enumerate(filepaths):
        # Read in Data and generate SpectrumModel
        Model = SpectrumModel.ReadInModelData(path, ZArray[count])
        # Convolve and Resample
        Model.ConvolveModel(FWHM=FWHM, PixelRes=Res)
        Model.ResampleFluxToGrid(ModelWavelengths)
        # Add Model to Array
        ModelArray = np.append(ModelArray, Model)

    LowerModel = ModelArray[1]
    UpperModel = ModelArray[2]
    Z_Int = np.power(10,0.5*(np.log10(LowerModel.Metallicity) + np.log10(UpperModel.Metallicity)))
    InterpolatedModel = SpectrumModel.InterpolateModels(Z_Int,ModelArray, Linear=True)

    fig, ax = plt.subplots(figsize=(8,8))
    NewModelArray = np.array([LowerModel, InterpolatedModel, UpperModel])
    Metallicities = [LowerModel.Metallicity, InterpolatedModel.Metallicity, UpperModel.Metallicity]
    Colors = ['#D62839', '#00916E', '#FF9F1C']
    for count, Model in enumerate(NewModelArray):
        ax.plot(Model.Wavelengths, Model.Fluxes, lw=.5, color=Colors[count], label='$\log Z_{*}$ = %.3f' % np.log10(Metallicities[count]))
    ax.legend(loc='upper right')
    ax.set_xlim(1250, 2000)
    ax.set_ylim(0.2e7, 0.9e7)
    ax.set_xlabel('Wavelength / Å')
    ax.set_ylabel('Model Flux')
    plt.savefig('FINPLOT/InterpolatedModels.png')
