# Relevant Imports
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib import rc
import emcee
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)

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

    Root = 'MassStacks'

    mode = 1 if Root == 'MassStacks' else 2
    #mode = 3

    for model in ['S', 'B']:
        if model == 'B':
            break

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
            MassData = Table.read(f'{Root}/Results/MassData/m{stack}_MassData.dat', format='ascii')
            StackMassDist = MassData['logMass'] - 10

            SFRData = Table.read(f'{Root}/Results/SFRData/m{stack}_SFRData.dat', format='ascii')
            StackSFRDist = SFRData['logSFR']

            # Determine Metallicity Distribution
            with open(f'{Root}/Results/{ModelType}/m{stack}Results.txt', 'r') as file:
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
        logMErrors = (MassArray[0], MassArray[2])
        logZErrors = (MetallicityArray[0], MetallicityArray[2])
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
        flat_samples = sampler.get_chain(discard=int(tau[0]*2), thin=int(tau[0]/2), flat=True)
        m_dist = np.percentile(flat_samples[:,0],[16,50,84])
        m_diff = np.diff(m_dist)
        print(m_dist, m_diff)
        c_dist = np.percentile(flat_samples[:,1],[16,50,84])
        c_diff = np.diff(c_dist)
        print(c_dist, c_diff)
        lowerlimits = [False]*7
        upperlimits = [True] + [False]*6

        m, c = m_dist[1], c_dist[1]

        C19m = 0.27
        C19c = -0.68

        fig, ax = plt.subplots(figsize=(8,6))
        ax.errorbar(logM, logZ, xerr=logMErrors, yerr=logZErrors, fmt='x', capsize=2, lw=.5, color='#D62839', label='FF Mass Stacks', uplims=upperlimits, lolims=lowerlimits)
        xlims = ax.get_xlim()
        xvals = np.linspace(xlims[0], xlims[1], 1000)
        ax.plot(xvals, m*xvals + c, lw=.5, ls='dashed', color='#D62839', label='FF MZR')
        ax.plot(xvals, C19m*xvals+C19c, lw=.5, ls='dashed', color='#90137E', label='C21 MZR')

        # Root
        Root = 'ContinuumFit_MassStacks'
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
            StackMassDist = MassData['logMass'] - 10

            SFRData = Table.read(f'{Root}/Results/SFRData/m{stack}_SFRData.dat', format='ascii')
            StackSFRDist = SFRData['logSFR']

            # Determine Metallicity Distribution
            with open(f'{Root}/Results/{ModelType}/m{stack}Results.txt', 'r') as file:
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
        logMErrors = (MassArray[0], MassArray[2])
        logZErrors = (MetallicityArray[0], MetallicityArray[2])
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
        flat_samples = sampler.get_chain(discard=int(tau[0]*2), thin=int(tau[0]/2), flat=True)
        m_dist = np.percentile(flat_samples[:,0],[16,50,84])
        m_diff = np.diff(m_dist)
        print(m_dist, m_diff)
        c_dist = np.percentile(flat_samples[:,1],[16,50,84])
        c_diff = np.diff(c_dist)
        print(c_dist, c_diff)
        lowerlimits = [False]*7
        upperlimits = [True] + [False]*6

        m, c = m_dist[1], c_dist[1]
        ax.errorbar(logM, logZ, xerr=logMErrors, yerr=logZErrors, fmt='x', capsize=2, lw=.5, color='#00916E', label='CN Mass Stacks', uplims=upperlimits, lolims=lowerlimits)
        ax.plot(xvals, m*xvals + c, lw=.5, ls='dashed', color='#00916E', label='CN MZR')

        ax.set_xlabel('$\log (M/M_{\odot})$')
        ax.set_ylabel('$\log (Z/Z_{\odot})$')
        ax.legend(loc='lower right')
        plt.tight_layout()
        plt.savefig(f'FINPLOT/FinalFigure.png')
