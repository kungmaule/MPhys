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
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(8,8), gridspec_kw={'height_ratios': [3,3,2]})

    #PertVar = np.zeros(0)
    #BinVar = np.zeros(0)

    # Iterate over all Stacks
    for count, stack in enumerate(FilePaths):
        # Fit Spline
        StackModel = SpectrumModel.ReadInTestData(stack)
        Errors = StackModel.Errors
        Fluxes = StackModel.Fluxes
        StackSpline = StackModel.DetermineContinuumSpline()
        Wl, SplFlux = StackSpline.Wavelengths, StackSpline.SplineFluxes
        PertVar = np.copy(StackSpline.SplineFluxes)
        break

    ax1.plot(StackSpline.Wavelengths, StackSpline.Fluxes * 1e19, color='black', lw=.5, label='Data')
    # Perturbed
    for i in range(499):
        PertModel = SpectrumModel.ReadInTestData(stack)
        x = np.random.normal(loc=0, scale=Errors)
        PertModel.Fluxes = PertModel.Fluxes + x
        PertSpline = PertModel.DetermineContinuumSpline()
        ax1.plot(PertSpline.Wavelengths, PertSpline.SplineFluxes * 1e19, color='#3C6997', alpha=0.125, lw=.5)
        PertVar = np.vstack((PertVar, PertSpline.SplineFluxes))

    #Plot all mass stack splines
    FileDirectory = f'./{Directory}/Stacks/{StackID}'
    PlotPath = f'./{Directory}/Plots'
    FilePaths = [os.path.join(FileDirectory, name) for name in os.listdir(FileDirectory)]
    for count, stack in enumerate(FilePaths):
        # Fit Spline
        StackModel = SpectrumModel.ReadInTestData(stack)
        StackSpline = StackModel.DetermineContinuumSpline()

        # Plot spectrum if unperturbed
        if count == 0:
            ax2.plot(StackSpline.Wavelengths, StackSpline.Fluxes * 1e19, color='black', lw=.5, label='Data')
            Wl, SplFlux = StackSpline.Wavelengths, StackSpline.SplineFluxes
            BinVar = np.copy(StackSpline.Fluxes)
        else:
            # Plot Spline
            ax2.plot(StackSpline.Wavelengths, StackSpline.SplineFluxes * 1e19, color='#3C6997', alpha=0.125, lw=.5)
            BinVar = np.vstack((BinVar, StackSpline.SplineFluxes))
        # Update Output
        print(f'[ Stack {count+1}/{len(FilePaths)} fitted. Overall Completion: {100*(count+1)/len(FilePaths)} % ]', end='\r')
    ax2.plot(Wl, SplFlux, color='#D62839', lw=.5, alpha=1, label='Unperturbed Spline')

    plt.xlim(StackSpline.Wavelengths[0], StackSpline.Wavelengths[-1])
    ax1.plot(Wl, SplFlux * 1e19, color='#D62839', lw=.5, alpha=1, label='Unperturbed Spline')
    WindowData = Table.read('Parameter Files/Windows/RixWindows.txt', format='ascii')
    WindowsStart = WindowData['wlmin']
    WindowsStop = WindowData['wlmax']
    for i in range(WindowsStart.shape[0]):
        ax1.axvspan(WindowsStart[i], WindowsStop[i], color='#C0C0C0', alpha=0.2)
        ax2.axvspan(WindowsStart[i], WindowsStop[i], color='#C0C0C0', alpha=0.2)
    ax3.set_xlabel('Wavelengths / Å')
    ax1.set_ylabel('Fluxes (Error Pert) $\\times 10^{-19}$')
    ax2.set_ylabel('Fluxes (Multi Stack) $\\times 10^{-19}$')
    ax2.set_ylim(ax1.get_ylim())
    ax1.legend(loc='lower left')
    ax2.legend(loc='lower left')
    plt.subplots_adjust(wspace=0, hspace=0)

    # ax 3
    Variance_Bins = np.empty(0)
    Variance_Pert = np.empty(0)
    print(BinVar.shape)
    for wl in range(BinVar.shape[1]):
        Variance_Bins = np.append(Variance_Bins, np.std(BinVar[:,wl]))
        Variance_Pert = np.append(Variance_Pert, np.std(PertVar[:,wl]))
    #Variance_Bins = np.var(BinVar, axis=0)
    #Variance_Pert = np.var(PertVar, axis=0)
    print(Variance_Bins)
    print(Variance_Pert)
    Bin_Mean = 1e20 * np.mean(Variance_Bins)
    Pert_Mean = 1e20 * np.mean(Variance_Pert)
    ax3.plot(Wl, Variance_Bins*1e20, color='black', lw=.5, label=f'Stack Mean: {Bin_Mean.round(3)}')
    ax3.plot(Wl, Variance_Pert*1e20, color='#D62839', lw=.5, label=f'Perturbation Mean: {Pert_Mean.round(3)}')
    ax3.set_ylabel('Spline Variance $\\times 10^{-20}$')
    ax3.legend(loc='upper left')
    plt.savefig(f'{PlotPath}/{StackID}_SplineVariations.png')
    # Plot Unperturbed Spectra for stack
    # Read in All stacks and plot continuum splines
