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
    if model.upper() == 'S':
        ModelType = 'S99'
        ZMin, ZMax = -3.0, -2.0
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

        StackData = Table.read(f'{Root}/Results/Distributions/m{stack}_disributions.dat', format='ascii')
        StackMassDist = StackData['lmass']
        StacklogZDist = StackData['lz']
        StackSFRDist = StackData['lsfr']

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

    # Simplify Error input
    logMErrors = (MassArray[0], MassArray[2])
    logZErrors = (MetallicityArray[0], MetallicityArray[2])
    logM = MassArray[1]
    logZ = MetallicityArray[1]

    # Initial Estimate
    grad, con = np.polyfit(logM, logZ, 1)

    MassArray = np.empty(shape=(3,7))
    MetallicityArray = np.empty(shape=(3,7))
    SFRArray = np.empty(shape=(3,7))

    ZErr, MErr = np.empty(7), np.empty(7)

    fig, ax = plt.subplots(figsize=(6,6))
    ax.errorbar(logM, logZ, xerr=logMErrors, yerr=logZErrors, fmt='x', lw=.5, capsize=2, color='#D62839', label='Full Fit', uplims=upperlimits, lolims=lowerlimits)
    FitX = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 1000)
    FitY = FitX*grad + con

    ax.plot(FitX, FitY, ls='--', color='#D62839', label='Fitted Linear Relation (Full Fit)', lw=.5)

    # CN

    for stack in range(1,8):

        StackData = Table.read(f'ContinuumFit_{Root}/Results/Distributions/m{stack}_disributions.dat', format='ascii')
        StackMassDist = StackData['lmass']
        StacklogZDist = StackData['lz']
        StackSFRDist = StackData['lsfr']

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

    # Simplify Error input
    logMErrors = (MassArray[0], MassArray[2])
    logZErrors = (MetallicityArray[0], MetallicityArray[2])
    logM = MassArray[1]
    logZ = MetallicityArray[1]

    # Initial Estimate
    grad, con = np.polyfit(logM, logZ, 1)


    # Plot Figure
    ax.errorbar(logM, logZ, xerr=logMErrors, yerr=logZErrors, fmt='x', lw=.5, capsize=2, color='#00916E', label='CN Fit', uplims=upperlimits, lolims=lowerlimits)
    # Fit
    FitX = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 1000)
    FitY = FitX*grad + con

    ax.plot(FitX, FitY, ls='--', color='#00916E', label='Fitted Linear Relation (CN)', lw=.5)
    #ax.fill_between(FitX, Fit_up, Fit_dw, color='blue', label='5-Sigma Interval', alpha=0.25)
    plt.ylabel('$\log (Z_{*}$)')
    plt.xlabel('$\log (M/M_{\odot}$)')
    plt.title('Plot of Mass-Metallicity Relationship')
    plt.legend(loc='upper left')
    plt.savefig(f'{Root}/Plots/{ModelType}/MZRvsCN_MZR.png')
    print('Gradient: ', grad, 'Constant:', con)
