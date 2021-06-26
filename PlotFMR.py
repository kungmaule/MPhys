# Relevant Imports
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib import rc
import dynesty as dyn
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)

##### MAIN #####################################################################

if __name__=="__main__":

    Root = sys.argv[1]
    model = sys.argv[2]
    mode = int(sys.argv[3])
    if model.upper() == 'S':
        ModelType = 'S99'
        #ZMin, ZMax = -3.0, -2.0
        lowerlimits = [False]*7
        upperlimits = [True] + [False]*6
    else:
        ModelType = 'BPASS'
        ZMin, ZMax = -4.0, -2.2
        lowerlimits = [False]*7
        upperlimits = [False]*7
    MassArray = np.empty(shape=(3,7,2))
    MetallicityArray = np.empty(shape=(3,7,2))
    SFRArray = np.empty(shape=(3,7,2))

    ZErr, MErr = np.empty(shape=(7,2)), np.empty(shape=(7,2))


    for stack in range(1,8):

        for sfr in range(1,3):
            StackMassDist = np.empty(0)
            StacklogZDist = np.empty(0)
            StackSFRDist = np.empty(0)

            # Determine Mass Distribution
            MassData = Table.read(f'{Root}/Results/MassData/m{stack}_sfr{sfr}_MassData.dat', format='ascii')
            StackMassDist = MassData['logMass']

            SFRData = Table.read(f'{Root}/Results/SFRData/m{stack}_sfr{sfr}_SFRData.dat', format='ascii')
            StackSFRDist = SFRData['logSFR']

            # Determine Metallicity Distribution
            with open(f'{Root}/Results/{ModelType}/m{stack}sfr{sfr}Results.txt', 'r') as file:
                for count, line in enumerate(file):
                    Tokens = line.split(', ')
                    StacklogZDist = np.append(StacklogZDist, np.log10(float(Tokens[2])/0.0142))

            # Save + Plot Distributions
            OutputTable = Table([StackMassDist, StackSFRDist, StacklogZDist], names=['lmass', 'lsfr', 'lz'])
            OutputTable.write(f'{Root}/Results/{ModelType}Distributions/m{stack}sfr{sfr}_disributions.dat', format='ascii', overwrite=True)


            # Take Means and Errors on Distribution
            logM = np.percentile(StackMassDist, [16, 50, 84])
            logMDifference = np.diff(logM)
            logSFR = np.percentile(StackSFRDist, [16, 50, 84])
            logSFRDifference = np.diff(logSFR)
            if stack == 1 and ModelType == 'S99':
                val = np.percentile(StacklogZDist, [68])
                err = np.percentile(StacklogZDist, [50])
                logZ = np.array([err,val,err])
            else:
                logZ = np.percentile(StacklogZDist, [16, 50, 84])
                logZDifference = np.diff(logZ)

            # Record Values
            MassArray[0, stack-1, sfr-1] = logM[1] - logM[0]
            MassArray[1, stack-1, sfr-1] = logM[1]
            MassArray[2, stack-1, sfr-1] = logM[2] - logM[1]
            MetallicityArray[0, stack-1, sfr-1] = logZ[1] - logZ[0]
            MetallicityArray[1, stack-1, sfr-1] = logZ[1]
            MetallicityArray[2, stack-1, sfr-1] = logZ[2] - logZ[1]
            SFRArray[0, stack-1, sfr-1] = logSFR[1] - logSFR[0]
            SFRArray[1, stack-1, sfr-1] = logSFR[1]
            SFRArray[2, stack-1, sfr-1] = logSFR[2] - logSFR[1]

            if mode == 2:
                fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6,12))
                ax1.hist(StackMassDist, bins=50, color='black', label='logM')
                ax1.axvline(MassArray[1,stack-1, sfr-1], color='red', label=f'Median: {logM[1]}')
                ax2.hist(StacklogZDist, bins=50, color='black', label='logZ')
                ax2.axvline(MetallicityArray[1,stack-1, sfr-1], color='red', label=f'Median: {logZ[1]}')
                ax3.hist(StackSFRDist, bins=50, color='black', label='logZ')
                ax3.axvline(SFRArray[1,stack-1, sfr-1], color='red', label=f'Median: {logSFR[1]}')
                ax1.set_ylabel('Frequency')
                ax1.set_xlabel('Log (M/$M_{\odot}$)')
                ax2.set_ylabel('Frequency')
                ax2.set_xlabel('Log (Z/$Z_{\odot}$)')
                ax3.set_xlabel('Log (SFR)')
                ax3.set_ylabel('Frequency')
                ax1.set_title(f'Distributions for Stack {stack}')
                #ax2.set_xlim(ZMin, ZMax)
                plt.savefig(f'{Root}/Plots/{ModelType}/m{stack}sfr{sfr}_Distributions.png')

    if mode == 2:
        quit()
    # Simplify Error input
    LOWlogMErrors = (MassArray[0,:,0], MassArray[2,:,0])
    LOWlogZErrors = (MetallicityArray[0,:,0], MetallicityArray[2,:,0])
    LOWlogM = MassArray[1,:,0]
    LOWlogZ = MetallicityArray[1,:,0]

    HIGHlogMErrors = (MassArray[0,:,1], MassArray[2,:,1])
    HIGHlogZErrors = (MetallicityArray[0,:,1], MetallicityArray[2,:,1])
    HIGHlogM = MassArray[1,:,1]
    HIGHlogZ = MetallicityArray[1,:,1]

    # Initial Estimate
    lowgrad, lowcon = np.polyfit(LOWlogM, LOWlogZ, 1)
    highgrad, highcon = np.polyfit(HIGHlogM, HIGHlogZ, 1)

    # Plot Figure
    fig, ax = plt.subplots(figsize=(6,6))
    ax.errorbar(LOWlogM, LOWlogZ, xerr=LOWlogMErrors, yerr=LOWlogZErrors, fmt='none', lw=.5, capsize=2, color='black', label='LOW SFR', uplims=upperlimits, lowlims=lowerlimits)
    ax.errorbar(HIGHlogM, HIGHlogZ, xerr=HIGHlogMErrors, yerr=HIGHlogZErrors, fmt='none', lw=.5, capsize=2, color='red', label='HIGH SFR', uplims=upperlimits, lowlims=lowerlimits)
    # Fit
    FitX = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 1000)
    LOWFitY = FitX*lowgrad + lowcon
    HIGHFitY = FitX*highgrad + highcon

    ax.plot(FitX, LOWFitY, ls='--', color='black', label='LOW SFR Fit', lw=.5)
    ax.plot(FitX, HIGHFitY, ls='--', color='red', label='HIGH SFR Fit', lw=.5)
    #ax.fill_between(FitX, Fit_up, Fit_dw, color='blue', label='5-Sigma Interval', alpha=0.25)
    plt.ylabel('Log (Z/$Z_{\odot}$)')
    plt.xlabel('Log (M/$M_{\odot}$)')
    plt.legend(loc='best')
    plt.savefig(f'{Root}/Plots/{ModelType}/SimpleFitFMR.png')
    print('Gradient: ', lowgrad, 'Constant:', lowcon)
    print('Gradient: ', highgrad, 'Constant:', highcon)
    plt.show()
