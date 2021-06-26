# Relevant Imports
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib import rc
import emcee
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=11)

##### METHODS ##################################################################

def log_likelihood(theta, x, y, yerr):
    m, c = theta
    model = m * x + c
    sigma2 = yerr**2
    return -0.5 * np.sum((y-model)**2 / sigma2)

def log_prior(theta):
    m, c = theta
    if 0 <= m <= 1:
        if -10 <= c <= 0:
            return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    if not np.isfinite(log_prior(theta)):
        return -np.inf
    return log_prior(theta) + log_likelihood(theta, x, y, yerr)

##### MAIN #####################################################################

if __name__=="__main__":

    Root = sys.argv[1]
    mode = int(sys.argv[2])
    #mode = 3
    fig, ax = plt.subplots(2, 2, figsize=(10,8), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    for model in ['S', 'B']:
        for sfr in [1,2]:

            if model.upper() == 'S':
                ModelType = 'S99'
                if mode == 1:
                    ZMin, ZMax = -1.2, -0.4
                else:
                    ZMin, ZMax = -1.2, 0.3
                lowerlimits = [False]*7
                upperlimits = [True] + [False]*6
                col = '#D62839'
            else:
                ModelType = 'BPASS'
                if mode == 1:
                    ZMin, ZMax = -2.0, -0.7
                else:
                    ZMin, ZMax = -3.1, -0.4
                lowerlimits = [False]*7
                upperlimits = [False]*7
                col = '#00916E'
            if mode == 3:
                if model == 'B':
                    break
                else:
                    #ModelType = 'Simulated (S99)'
                    ID = 'Simulated S99'
                    col = '#D62839'

            MassArray = np.empty(shape=(3,7))
            MetallicityArray = np.empty(shape=(3,7))
            SFRArray = np.empty(shape=(3,7))

            ZErr, MErr = np.empty(7), np.empty(7)

            for stack in range(1,8):

                StackMassDist = np.empty(0)
                StacklogZDist = np.empty(0)
                StackSFRDist = np.empty(0)

                # Determine Mass Distribution
                MassData = Table.read(f'{Root}/Results/MassData/m{stack}_sfr{sfr}_MassData.dat', format='ascii')
                StackMassDist = MassData['logMass']
                if mode == 2:
                    StackMassDist -= 10

                SFRData = Table.read(f'{Root}/Results/SFRData/m{stack}_sfr{sfr}_SFRData.dat', format='ascii')
                StackSFRDist = SFRData['logSFR']

                # Determine Metallicity Distribution
                with open(f'{Root}/Results/{ModelType}/m{stack}sfr{sfr}Results.txt', 'r') as file:
                    for count, line in enumerate(file):
                        Tokens = line.split(', ')
                        StacklogZDist = np.append(StacklogZDist, np.log10(float(Tokens[2])/0.0142))

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

            # Simplify Error input
            logMErrors = [MassArray[0], MassArray[2]]
            logZErrors = np.array([MetallicityArray[0], MetallicityArray[2]])
            logM = MassArray[1]
            logZ = MetallicityArray[1]

            XFittingErrors = np.maximum(MassArray[0], MassArray[2])
            YFittingErrors = np.maximum(MetallicityArray[0], MetallicityArray[2])

            # Initial Estimate
            m_guess, c_guess  = np.polyfit(logM, logZ, 1)

            pos = np.array([m_guess, c_guess]) + 1e-4*np.random.randn(16,2)
            nwalkers, ndim = pos.shape

            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(logM, logZ, YFittingErrors))
            sampler.run_mcmc(pos, 10000, progress=True)
            tau = sampler.get_autocorr_time()
            print(tau)
            #discard = int(input('Discard: '))
            #thin = int(input('Thin: '))
            discard, thin = int(tau[0]*2), int(tau[0]/2)
            flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=True)
            m_dist = np.percentile(flat_samples[:,0],[16,50,84])
            m_diff = np.diff(m_dist)
            print(ModelType)
            print(sfr, m_dist, m_diff)
            c_dist = np.percentile(flat_samples[:,1],[16,50,84])
            c_diff = np.diff(c_dist)
            print(sfr, c_dist, c_diff)

            m, c = m_dist[1], c_dist[1]
            if sfr==1:
                sfrlabel = 'Low'
                color='#00916E'
            else:
                sfrlabel = 'High'
                color='#D62839'

            if mode == 2 and sfr==2 and ModelType=='B':
                quit()
            if model == 'S':
                ax[0,0].errorbar(logM, logZ, xerr=logMErrors, yerr=logZErrors, fmt='none', lw=.5, capsize=2, color=color, label=f'{sfrlabel} SFR', uplims=upperlimits, lowlims=lowerlimits)
                FitX = np.linspace(ax[0,0].get_xlim()[0], ax[0,0].get_xlim()[1], 1000)
                UpperYlim = ax[0,0].get_ylim()[1]
                ax[0,0].plot(FitX, m_dist[1]*FitX + c_dist[1], ls='--', color=color, label=f'{sfrlabel} SFR Fit', lw=.5)
                if sfr == 1:
                    a, b = 0.29, -3.57
                    ax[0,0].plot(FitX, a*FitX+b, linestyle='dotted', lw=.5, color='black', label='MZR Fit' )
                    Ssfr1 = logZ
                    Ssfr1err = logZErrors
                else:
                    Error = np.sqrt(logZErrors**2 + Ssfr1err**2)
                    ax[1,0].errorbar(logM, logZ - Ssfr1, yerr=Error, fmt='x', lw=.5, capsize=2, color='#90137E')
                    ax[1,0].set_ylabel('$\Delta \log (Z/Z_{\odot})$')
                    ax[1,0].axhline(0.0, color='black', lw=.5)
            if model=='B':
                ax[0,1].errorbar(logM, logZ, xerr=logMErrors, yerr=logZErrors, fmt='none', lw=.5, capsize=2, color=color, label=f'{sfrlabel} SFR', uplims=upperlimits, lowlims=lowerlimits)
                FitX = np.linspace(ax[0,1].get_xlim()[0], ax[0,1].get_xlim()[1], 1000)
                ax[0,1].plot(FitX, m_dist[1]*FitX + c_dist[1], ls='--', color=color, label=f'{sfrlabel} SFR Fit', lw=.5)
                LowerYlim = ax[0,0].get_ylim()[0]
                if sfr == 1:
                    a, b = 0.32, -4.22
                    ax[0,1].plot(FitX, a*FitX+b, linestyle='dotted', lw=.5, color='black', label='MZR Fit' )
                    Ssfr1 = logZ
                    Ssfr1err = logZErrors
                else:
                    Error = np.sqrt(logZErrors**2 + Ssfr1err**2)
                    ax[1,1].errorbar(logM, logZ - Ssfr1, yerr=Error, fmt='x', lw=.5, capsize=2, color='#90137E')
                    ax[1,1].set_ylabel('$\Delta \log (Z/Z_{\odot})$')
                    ax[1,1].axhline(0.0, color='black', lw=.5)

    #ax[0,0].set_title('S99 Fits')
    ax[0,1].set_ylabel('Log (Z/$Z_{\odot}$)')
    #ax[0,0].set_ylim(-1.1,-0.6)
    #ax[0,1].set_ylim(-2.0,-0.8)
    ax[0,0].set_ylabel('$\log (Z/Z_{\odot})$')
    ax[1,0].set_xlabel('$\log (M/M_{\odot})$')
    ax[1,1].set_xlabel('$\log (M/M_{\odot})$')
    #ax[0,1].set_ylim(ax[0,0].get_ylim())
    #ax[1,1].set_ylim(ax[1,0].get_ylim())
    plt.subplots_adjust(hspace=0)
    #ax[0,1].set_title('BPASS Fits')
    #ax[0,0].set_title('Fundamental Metallicity Relationship')
    ax[0,0].legend(loc='best', title='$\\textbf{S99}$', fontsize=10)
    ax[0,1].legend(loc='best', title='$\\textbf{BPASS}$', fontsize=10)
    plt.savefig(f'FINPLOT/FMR_Plot.png')
