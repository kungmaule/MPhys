# Imports from libraries
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import rc
rc('text', usetex=True)
import os
import sys
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
# Imports from Classes
from SpectrumModel import SpectrumModel, Spline, GenerateColours

##### MAIN METHOD ##############################################################

if __name__=="__main__":

    Root = sys.argv[1]
    Mod = sys.argv[2]

    # Define Necessary Constants
    ModelSplineWindows = 'Parameter Files/Windows/RixWindows(Max2000).txt'
    DataSplineWindows = 'Parameter Files/Windows/RixWindows(Max2000).txt'
    ChiSquareWindows = 'Parameter Files/Windows/SteidelWindows.txt'

    # Generate Model Array
    if Mod == 'S':
        ZArray = np.array([0.001, 0.002, 0.008, 0.014, 0.040])
        FWHM, Res = 3.0, 0.4 # 0.4Å Resolution
        folderpath = r"./S99"
        ResultPath = 'S99'
    elif Mod == 'B':
        ZArray = np.array([0.00001, 0.0001, 0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.010, 0.014, 0.020, 0.030])
        FWHM, Res = 3.0, 1.0 # 1Å Resolution
        folderpath = r"./BPASS"
        ResultPath = 'BPASS'
    else:
        print("Invalid Model input. <S> for S99, <B> for BPASS.")
        quit()
    #ZArraySolar = (ZArray/0.0142).round(2)
    filepaths = [os.path.join(folderpath, name) for name in os.listdir(folderpath)]
    ModelArray = np.empty(0)
    Stepper = 0
    ModelWavelengths = np.arange(1000, 2000, 1)
    for count, path in enumerate(filepaths):
        # Read in Data and generate SpectrumModel
        Model = SpectrumModel.ReadInModelData(path, ZArray[count])
        # Convolve and Resample
        Model.ConvolveModel(FWHM=FWHM, PixelRes=Res)
        Model.ResampleFluxToGrid(ModelWavelengths)
        # Add Model to Array
        ModelArray = np.append(ModelArray, Model)

    #Trues = [0.092, 0.113, 0.134, 0.156, 0.178, 0.199, 0.220, 0.242, 0.263, 0.285, 0.306, 0.328, 0.350, 0.371, 0.392, 0.414, 0.436, 0.457, 0.478, 0.5]
    #Trues = [0.113, 0.156, 0.199, 0.242, 0.285, 0.328, 0.371, 0.414, 0.457, 0.500]
    #Trues = np.array([-2.98, -2.78, -2.78, -2.70, -2.71, -2.56, -2.42])
    #Trues = np.power(10, Trues)

    DataFile = open(f'{Root}/Data/Results.txt')
    for line in DataFile:
        Tokens = line.split(', ')
        TestIndex = int(Tokens[0])
        OutZ = float(Tokens[2])
        #TrueZ = Trues[TestIndex-1]

        # Define ZErr
        LowZ = float(Tokens[1])
        HighZ = float(Tokens[3])

        # Read in Test Data
        SpectrumPath = f"{Root}/Spectra/vandels-m{TestIndex}-p0-allz.dat"
        # Generate Model based upon test data
        StackModel = SpectrumModel.ReadInTestData(SpectrumPath)
        # Fit Spline to Model
        StackSpline = StackModel.DetermineContinuumSpline()

        # Colours
        CArr = GenerateColours()

        # Generate Plot
        fig = plt.figure(figsize=(14, 7))
        gs = fig.add_gridspec(nrows=3,ncols=2, width_ratios=[3,2], height_ratios=[5,5,2])
        axSpl = fig.add_subplot(gs[0, :-1])
        axBFM = fig.add_subplot(gs[1, :-1], sharex=axSpl)
        axErr = fig.add_subplot(gs[-1, :-1], sharex=axSpl)
        axPos = fig.add_subplot(gs[:, -1])
        plt.subplots_adjust(wspace=0.1, hspace=0.0)
        plt.setp(axSpl.get_xticklabels(), visible=False)
        plt.setp(axBFM.get_xticklabels(), visible=False)

        # Plot Posterior Data
        Pos = Table.read(f'{Root}/Data/m{TestIndex}_Posterior.dat', format='ascii')
        PosData = Pos['pos']
        PosMean = np.mean(PosData)
        axPos.hist(PosData, 50, density=False, color=CArr[0])
        axPos.axvline(PosMean, color=CArr[6], label='Mean: %.3f' % PosMean, lw=.5)
        #axPos.axvline(TrueZ, color=CArr[4], label='True: %.3f' % TrueZ, lw=.5)
        axPos.set_xlabel('Metallicity')
        axPos.set_ylabel('Posterior Probability')
        axPos.set_title(f'Posterior Probability of NS: m{TestIndex}')
        axPos.legend(loc='upper right')

        # Plot Spline Data
        axSpl.plot(StackSpline.Wavelengths, StackSpline.Fluxes, color=CArr[1], lw=.4, label='Spectrum Data')
        axSpl.plot(StackSpline.Wavelengths, StackSpline.SplineFluxes, color=CArr[2], lw=.4, label='Spline')
        axSpl.set_xlim(StackSpline.Wavelengths[0], StackSpline.Wavelengths[-1])
        axSpl.legend(loc='upper right')
        axSpl.set_ylabel('Flux / Spline')
        axSpl.set_title('Spline, CN Data vs Best Fit Model, Continuum Normalised Error')
        WindowData = Table.read(DataSplineWindows, format='ascii')
        WindowsStart = WindowData['wlmin']
        WindowsStop = WindowData['wlmax']
        for i in range(WindowsStart.shape[0]):
            axSpl.axvspan(WindowsStart[i], WindowsStop[i], color='#C0C0C0', alpha=0.2)


        # Plot BFM Data
        ZArr = np.array([LowZ, OutZ, HighZ])
        ZMods = np.empty(shape=ZArr.shape[0], dtype=object)
        for i in range(ZArr.shape[0]):
            Z = ZArr[i]
            if Z in ZArray:
                Model = ModelArray[np.where(ZArray==OutZ)[0][0]]
            else:
                Model = SpectrumModel.InterpolateModels(OutZ, ModelArray, Linear=True)
            # Convolve to correct Resolution
            Model.ConvolveModel(FWHM=3.0, PixelRes=0.4)
            # Determine Relevant Continuum Spline
            SplineModel = Model.DetermineContinuumSpline()
            # Resample onto Wavelength Grid
            SplineModel.ResampleRelativeFluxToGrid(StackSpline.Wavelengths)
            ZMods[i] = SplineModel
        axBFM.plot(StackSpline.Wavelengths, StackSpline.RelativeFluxes, color=CArr[1], lw=.4, label='CN Flux')
        axBFM.plot(ZMods[1].Wavelengths, ZMods[1].RelativeFluxes, color=CArr[2], lw=.4, label='Best Fit Model: %.3f $Z\\textsubscript{\(\odot\)}$' % ZMods[1].Metallicity)
        MinFlux = np.minimum(ZMods[2].RelativeFluxes, ZMods[0].RelativeFluxes)
        MaxFlux = np.maximum(ZMods[2].RelativeFluxes, ZMods[0].RelativeFluxes)
        axBFM.fill_between(StackSpline.Wavelengths, MinFlux, MaxFlux, alpha=0.7, color=CArr[2], label='Best Fit Model Error')
        axBFM.hlines(1.0, StackSpline.Wavelengths[0], StackSpline.Wavelengths[-1], color='black', lw=.1)
        axBFM.set_ylabel('Continuum Normalised Flux')
        WindowData = Table.read(ChiSquareWindows, format='ascii')
        WindowsStart = WindowData['wlmin']
        WindowsStop = WindowData['wlmax']
        for i in range(WindowsStart.shape[0]-1):
            if i == 0:
                axBFM.axvspan(StackSpline.Wavelengths[0], WindowsStart[i], color='#C0C0C0', alpha=0.2)
                axErr.axvspan(StackSpline.Wavelengths[0], WindowsStart[i], color='#C0C0C0', alpha=0.2)
            elif i == WindowsStart.shape[0]-1:
                axBFM.axvspan(WindowsStop[i], StackSpline.Wavelengths[-1], color='#C0C0C0', alpha=0.2)
                axErr.axvspan(WindowsStop[i], StackSpline.Wavelengths[-1], color='#C0C0C0', alpha=0.2)
            else:
                axBFM.axvspan(WindowsStop[i], WindowsStart[i+1], color='#C0C0C0', alpha=0.2)
                axErr.axvspan(WindowsStop[i], WindowsStart[i+1], color='#C0C0C0', alpha=0.2)


        axBFM.legend(loc='upper right')


        # Plot Error Data
        axErr.plot(StackSpline.Wavelengths, StackSpline.Errors, color=CArr[5], lw=.4, label='Spectral Error')
        axErr.set_ylabel('CN Error')
        axErr.set_xlabel('Wavelength / Å')
        axErr.legend(loc='upper right')


        plt.savefig(f'{Root}/Plots/SplBfmPos/m{TestIndex}.png')
        #axBFM, axErr = BFM(Data, Z, ZErr, ModelArray, Windows=ChiSquareWindows, TestIndex, TrueZ)

        #axSpline = ContinuumSpline(Data, Windows=DataSplineWindows)
        print(f'Test Spectrum {TestIndex} plotted.', end='\r')
        # END LOOP
    print("All Test Spectra successfully plotted.")
