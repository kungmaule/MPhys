################################################################################

"""             Thomas Stanton - Metallicity Fit NS Methodology            """

################################################################################

# Imports from libraries
from astropy.table import Table
import numpy as np
import os, sys

# Import Fitting libraries
import dynesty
from multiprocessing import Pool

# Imports from Classes
from SpectrumModel import SpectrumModel, Spline

##### MCMC METHODS #############################################################

def PriorTransform(u, ZMin, ZMax):
    return ZMin + u*(ZMax-ZMin)

def LogLikelihood(Z, Fluxes, Wavelengths, Errors, ModelArray, Mask):

    # Get Model Fluxes at Z
    ModelFluxes = InterpolateZModel(Z, ModelArray, Wavelengths)

    # Mask All Arrays
    Wavelengths = Wavelengths[Mask]
    MaskFlux = Fluxes[Mask]
    ModelMaskFluxes = ModelFluxes[Mask]
    MaskErrors = Errors[Mask]

    # Determins s^2
    sigma2 = MaskErrors**2

    # Return Log Likelihood
    return -0.5 * np.sum(((MaskFlux - ModelMaskFluxes)**2 / sigma2) + np.log(2*np.pi*sigma2))

def InterpolateZModel(Z, ModelArray, Wavelengths):

    # Define initial Metallicity Array
    ZArray = np.array([Model.Metallicity for Model in ModelArray])
    # Generate Model based on Metallicity
    Model = ModelArray[np.where(ZArray==Z)][0][0] if Z in ZArray else SpectrumModel.InterpolateModels(Z,ModelArray, Linear=True)
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

def FitMetallicity(stack_path, stack_id, pert_id, model_type, ResultsPath, Verbose=True): #m_id, p_id, model_type, Verbose=True):

    # Generate Model Array
    if model_type == 'S':
        ZArray = np.array([0.001, 0.002, 0.008, 0.014, 0.040])
        FWHM, Res = 3.0, 0.4 # 0.4Å Resolution
        folderpath = r"./S99"
        ModelPath = 'S99'
    elif model_type == 'B':
        ZArray = np.array([0.00001, 0.0001, 0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.010, 0.014, 0.020, 0.030])
        FWHM, Res = 3.0, 1.0 # 1Å Resolution
        folderpath = r"./BPASS"
        ModelPath = 'BPASS'
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
    #SpectrumPath = f'Sample/Spectra/vandels-m{m_id}-p0-allz.dat'

    # Generate Model based upon test data + fit spline
    StackModel = SpectrumModel.ReadInTestData(stack_path)
    StackSpline = StackModel.DetermineContinuumSpline()

    # Assign Values to easier to reference variables
    Flux = StackSpline.RelativeFluxes
    Wavelengths = StackSpline.Wavelengths
    Errors = StackSpline.Errors

    # Define Steidel Mask
    Mask = GenerateSteidelMask(Wavelengths)

    # Run Nested Sampling
    with Pool() as pool:
        sampler = dynesty.NestedSampler(LogLikelihood, PriorTransform, ndim=1, logl_args=(Flux, Wavelengths, Errors, ModelArray, Mask), ptform_args=(ZMin, ZMax), pool=pool, queue_size=6)
        sampler.run_nested(print_progress=Verbose)

    # Collect Results + Define Weights
    Results = sampler.results
    Weights = np.exp(Results['logwt'] - Results['logz'][-1])

    # Determine Metallicity and Posteriors
    metallicity = dynesty.utils.quantile(x=Results['samples'][:,0], q=[0.16, 0.50, 0.84], weights=Weights)
    samples_equal = dynesty.utils.resample_equal(Results['samples'], Weights)

    # Output if Verbose
    if Verbose == True: print(metallicity)

    # Record Metallicity Data
    with open(f'{ResultsPath}/{ModelPath}/{stack_id}Results.txt', 'a') as file:
        file.write(f'{pert_id}, {metallicity[0]}, {metallicity[1]}, {metallicity[2]}\n')
