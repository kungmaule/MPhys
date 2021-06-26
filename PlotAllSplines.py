# Numerical / System Imports
from astropy.table import Table
import numpy as np
import os, sys

# Graphing Imports
from matplotlib import colors, cm, gridspec, rc
import matplotlib.pyplot as plt
#from matplotlib import gridspec
from matplotlib import rc
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)
# Imports from Classes
from SpectrumModel import SpectrumModel, Spline

##### MAIN #####################################################################

if __name__=="__main__":

    if len(sys.argv) == 1:
        print(f'Usage: python {sys.argv[0]} <Directory> <Stack ID e.g. m1>')
        quit()

    # Read input
    Directory = sys.argv[1]
    StackID = sys.argv[2]

    # Read in Unperturbed Spectrum

    # Set up file pathing
    FileDirectory = f'./{Directory}/Stacks/{StackID}'
    PlotPath = f'./{Directory}/Plots'
    FilePaths = [os.path.join(FileDirectory, name) for name in os.listdir(FileDirectory)]

    # Set up Plot
    plt.figure(figsize=(12,6))

    # Iterate over all Stacks
    for count, stack in enumerate(FilePaths):
        # Fit Spline
        StackModel = SpectrumModel.ReadInTestData(stack)
        StackSpline = StackModel.DetermineContinuumSpline()

        # Plot spectrum if unperturbed
        if count == 0:
            plt.plot(StackSpline.Wavelengths, StackSpline.Fluxes, color='black', lw=.5, label='Data')
            Wl, SplFlux = StackSpline.Wavelengths, StackSpline.SplineFluxes
            plt.xlim(StackSpline.Wavelengths[0], StackSpline.Wavelengths[-1])
        else:
            # Plot Spline
            plt.plot(StackSpline.Wavelengths, StackSpline.SplineFluxes, color='#3C6997', alpha=0.125, lw=.5)
        # Update Output
        print(f'[ Stack {count+1}/{len(FilePaths)} fitted. Overall Completion: {100*(count+1)/len(FilePaths)} % ]', end='\r')
    plt.plot(Wl, SplFlux, color='#D62839', lw=.5, alpha=1, label='Unperturbed Spline')
    WindowData = Table.read('Parameter Files/Windows/RixWindows.txt', format='ascii')
    WindowsStart = WindowData['wlmin']
    WindowsStop = WindowData['wlmax']
    for i in range(WindowsStart.shape[0]):
        plt.axvspan(WindowsStart[i], WindowsStop[i], color='#C0C0C0', alpha=0.2)
    plt.xlabel('Wavelengths / Å')
    plt.ylabel('Fluxes')
    plt.title(f'Plot of Splines for {Directory}: {StackID}')
    plt.legend(loc='lower left')
    plt.savefig(f'{PlotPath}/{StackID}_AllSplines.png')
    # Plot Unperturbed Spectra for stack
    # Read in All stacks and plot continuum splines
