################################################################################

"""             Thomas Stanton - Metallicity Fit MCMC Methodology            """

################################################################################

# Imports from libraries
from astropy.table import Table
import numpy as np
import os, sys

# Import MCMC libraries
import emcee, corner
from multiprocessing import Pool

# Imports from Classes
from SpectrumModel import SpectrumModel, Spline

##### MCMC METHODS #############################################################

def LogLikelihood(Z, Fluxes, Wavelengths, Errors, ModelArray, Mask):

    # Get Model Fluxes at Z
    ModelFluxes = InterpolateZModel(Z, ModelArray, Wavelengths)
    Wavelengths = Wavelengths[Mask]
    MaskFlux = Fluxes[Mask]
    ModelMaskFluxes = ModelFluxes[Mask]
    MaskErrors = Errors[Mask]

    # Determins s^2
    sigma2 = MaskErrors**2

    # Return Log Likelihood
    return -0.5 * np.sum((MaskFlux - ModelMaskFluxes)**2 / sigma2 + np.log(2*np.pi*sigma2))

def LogPrior(Z, ZMin, ZMax):
    # Prior Function
    if ZMin <= Z <= ZMax:
        return 0.0
    return -np.inf

def LogProbability(Z, Fluxes, Wavelengths, Errors, ModelArray, Mask, ZMin, ZMax):
    lp = LogPrior(Z, ZMin, ZMax)
    if not np.isfinite(lp):
        return -np.inf
    return lp + LogLikelihood(Z, Fluxes, Wavelengths, Errors, ModelArray, Mask)

def InterpolateZModel(Z, ModelArray, Wavelengths):

    # Define initial Metallicity Array
    ZArray = np.array([Model.Metallicity for Model in ModelArray])
    # Generate Model based on Metallicity
    Model = ModelArray[np.where(ZArray==Z)][0]cd M   if Z in ZArray else SpectrumModel.InterpolateModels(Z,ModelArray, Linear=True)
    # Determine Continuum Spline
    SplineModel = Model.DetermineContinuumSpline()

    return SplineModel.RelativeFluxes

##### METHODS ##################################################################

def GenerateSteidelMask(Wavelengths):
    ContinuumWindows = Table.read('Parameter Files/Windows/SteidelWindows.txt', format='ascii')
    WindowWlMins = ContinuumWindows['wlmin']
    WindowWlMaxs = ContinuumWindows['wlmax']

    Mask = np.array([np.any((WindowWlMins <= wl) & (wl <= WindowWlMaxs)) for wl in Wavelengths])

    return Mask

##### MAIN METHOD ##############################################################

def FitMetallicity(rank_id, m_id, p_id, model_type):

    # Generate Model Array
    if model_type == 'S':
        ZArray = np.array([0.001, 0.002, 0.008, 0.014, 0.040])
        FWHM, Res = 3.0, 0.4 # 0.4Å Resolution
        folderpath = r"./S99"
        ResultPath = 'S99'
    elif model_type == 'B':
        ZArray = np.array([0.00001, 0.0001, 0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.010, 0.014, 0.020, 0.030])
        FWHM, Res = 3.0, 1.0 # 1Å Resolution
        folderpath = r"./BPASS"
        ResultPath = 'BPASS'
    ZMin, ZMax = ZArray[0], ZArray[-1]

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


    # Read in Spectrum Data
    #SpectrumPath = f'GALAXY_STACKS/MASS_STACKS/m{m_id}/vandels-m{m_id}-p{p_id}-allz.dat'
    SpectrumPath = f'MethodologyTests/Spectra/test_spectrum{m_id}.dat'
    # Generate Model based upon test data
    StackModel = SpectrumModel.ReadInTestData(SpectrumPath)
    # Fit Spline to Model
    StackSpline = StackModel.DetermineContinuumSpline()
    # Rescale Errors
    Flux = StackSpline.RelativeFluxes
    Wavelengths = StackSpline.Wavelengths
    Errors = StackSpline.Errors

    # Define Steidel Mask
    Mask = GenerateSteidelMask(Wavelengths)

    # Initialise Walker Data
    WalkerNumber = 6
    pos = np.empty(shape=(WalkerNumber, 1))
    positions = np.linspace(ZMin, ZMax, WalkerNumber)
    for count, position in enumerate(positions):
        pos[count] = position
    nwalkers, ndim = pos.shape

    # Initialise how many samples to run
    SampleNumber = 4000

    # Run Markov Chain
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, LogProbability, args=(Flux, Wavelengths, Errors, ModelArray, Mask, ZMin, ZMax), pool=pool)
        sampler.run_mcmc(pos, SampleNumber, progress=True)

    # Autocorrelation parameters
    #print(sampler.get_autocorr_time())
    tau = sampler.get_autocorr_time()
    print(tau)
    thin_tau = tau[0]/2
    disc_tau = tau[0]*2
    Thinned = int(thin_tau.round(0))
    Discarded = int(disc_tau.round(0))

    # Flatten Chain
    flat_samples = sampler.get_chain(discard=Discarded, thin=Thinned, flat=True)

    # Determine Mean and Errors
    mcmc = np.percentile(flat_samples[:,0], [16, 50, 84])
    difference = np.diff(mcmc)
    print(mcmc[1]/0.0142)

    # Record Data
    #print(f'Recording Data: {rank_id} {p_id}')
    with open(f'MethodologyTests/EmceeVsDynesty/EmceeResults.txt', 'a') as file:
        file.write(f'{int(m_id)-1}, {mcmc[1]}, {difference[0]}, {difference[1]} \n')


if __name__=="__main__":
    for m in range(20):
        rank = 100
        FitMetallicity(rank, m, 0, 'S')
