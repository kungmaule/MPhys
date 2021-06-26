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
    MassArray = np.empty(shape=(3,7))
    MetallicityArray = np.empty(shape=(3,7))
    SFRArray = np.empty(shape=(3,7))

    ZErr, MErr = np.empty(7), np.empty(7)

    for stack in range(1,8):

        StackMassDist = np.empty(0)
        StacklogZDist = np.empty(0)
        StackSFRDist = np.empty(0)

        # Determine Mass Distribution
        MassData = Table.read(f'{Root}/Results/MassData/m{stack}_MassData.dat', format='ascii')
        StackMassDist = MassData['logMass']

        SFRData = Table.read(f'{Root}/Results/SFRData/m{stack}_SFRData.dat', format='ascii')
        StackSFRDist = SFRData['logSFR']

        # Determine Metallicity Distribution
        with open(f'{Root}/Results/{ModelType}/m{stack}Results.txt', 'r') as file:
            for count, line in enumerate(file):
                Tokens = line.split(', ')
                StacklogZDist = np.append(StacklogZDist, np.log10(float(Tokens[2])/0.0142))

        # Save + Plot Distributions
        OutputTable = Table([StackMassDist, StackSFRDist, StacklogZDist], names=['lmass', 'lsfr', 'lz'])
        OutputTable.write(f'{Root}/Results/Distributions/{ModelType}/m{stack}_distributions.dat', format='ascii', overwrite=True)


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
        MassArray[0, stack-1] = logM[1] - logM[0]
        MassArray[1, stack-1] = logM[1]
        MassArray[2, stack-1] = logM[2] - logM[1]
        MetallicityArray[0, stack-1] = logZ[1] - logZ[0]
        MetallicityArray[1, stack-1] = logZ[1]
        MetallicityArray[2, stack-1] = logZ[2] - logZ[1]
        SFRArray[0, stack-1] = logSFR[1] - logSFR[0]
        SFRArray[1, stack-1] = logSFR[1]
        SFRArray[2, stack-1] = logSFR[2] - logSFR[1]

        ZErr[stack-1] = np.std(StacklogZDist)
        MErr[stack-1] = np.std(StackMassDist)

        if mode == 2:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,5))
            ax1.hist(StackMassDist, bins=50, color='black', label='logM')
            ax1.axvline(MassArray[1,stack-1], color='red', label=f'Median: {logM[1]}')
            ax2.hist(StacklogZDist, bins=50, color='black', label='logZ')
            ax2.axvline(MetallicityArray[1,stack-1], color='red', label=f'Median: {logZ[1]}')
            ax1.set_ylabel('Frequency')
            ax1.set_xlabel('Log (M/$M_{\odot}$)')
            ax2.set_ylabel('Frequency')
            ax2.set_xlabel('Log (Z/$Z_{\odot}$)')
            plt.savefig(f'{Root}/Plots/{ModelType}/m{stack}_Distributions.png')

    # Simplify Error input
    if mode == 2:
        quit()
    logMErrors = (MassArray[0], MassArray[2])
    logZErrors = (MetallicityArray[0], MetallicityArray[2])
    logM = MassArray[1]
    logZ = MetallicityArray[1]

    # Initial Estimate
    grad, con = np.polyfit(logM, logZ, 1)

    # Dynesty
    ndim = 2

    # Plot Figure
    fig, ax = plt.subplots(figsize=(6,6))
    ax.errorbar(logM, logZ, xerr=logMErrors, yerr=logZErrors, fmt='none', lw=.5, capsize=2, label='Data Points', uplims=upperlimits, lowlims=lowerlimits)

    # Fit
    FitX = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 1000)
    FitY = FitX*grad + con

    ax.plot(FitX, FitY, ls='--', color='black', label='Fitted Linear Relation')
    #ax.fill_between(FitX, Fit_up, Fit_dw, color='blue', label='5-Sigma Interval', alpha=0.25)
    plt.ylabel('Log (Z/$Z_{\odot}$)')
    plt.xlabel('Log (M/$M_{\odot}$)')
    plt.title('Plot of Mass-Metallicity Relationship')
    #plt.savefig(f'{Root}/Plots/{ModelType}/SimpleFitMZR.png')
    print('Gradient: ', grad, 'Constant:', con)
    plt.show()
