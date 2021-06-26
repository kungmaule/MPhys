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

def attenuation_salim_2018(wl, av, B, delta):
    """
    Returns A(lambda) for the Salim + 2018 dust attenuation
    prescription
    """

    x = 1.0 / wl
    rv_calz = 4.05
    k_calz = (2.659 * (-2.156 + 1.509*x - 0.198*x**2 + 0.011*x**3)) + rv_calz

    wl0 = 0.2175
    dwl = 0.035
    d_lam = (B * np.power(wl*dwl, 2)) / (np.power(wl**2-wl0**2, 2) + np.power(wl0*dwl, 2))

    rv_mod = rv_calz / ((rv_calz + 1)*(0.44/0.55)**delta - rv_calz)

    kmod = k_calz * (rv_mod/rv_calz) * (wl/0.55)**delta + d_lam

    return (kmod * av) / rv_mod

def PriorTransform(u, ZMin, ZMax):

    # Make copy of parameter set
    theta = np.copy(u)

    # Metallicity
    theta[0] = ZMin + u[0] * (ZMax - ZMin)
    # Av
    theta[1] = u[1]*5
    # B
    theta[2] = u[2]*5
    # Delta
    theta[3] = -1 + u[3]*2

    return theta

def LogLikelihood(theta, Fluxes, Wavelengths, Errors, ModelArray, Mask, ModelMask):

    # Unpack Parameters
    Z, Av, B, delta = theta

    # Get Model Fluxes at Z + Generate Salim Curve
    Z_Model = InterpolateZModel(Z, ModelArray)
    A_Lambda = attenuation_salim_2018(wl=Z_Model.Wavelengths/1.e4, av=Av, B=B, delta=delta)

    # Apply Dust
    ModelFluxes = Z_Model.Fluxes * np.power(10, -0.4*A_Lambda)

    # Mask All Arrays
    Wavelengths = Wavelengths[Mask]
    MaskFlux = Fluxes[Mask]
    ModelMaskFluxes = ModelFluxes[ModelMask]
    MaskErrors = Errors[Mask]
    # Determins s^2
    sigma2 = MaskErrors**2

    # Calculate Beta
    Beta = np.sum(ModelMaskFluxes*MaskFlux/sigma2) / np.sum(np.power(ModelMaskFluxes/MaskErrors,2))


    # Return Log Likelihood
    return -0.5 * np.sum(np.power(((Beta*ModelMaskFluxes)-MaskFlux)/MaskErrors,2) + np.log(2*np.pi*sigma2))

    #Beta*((MaskFlux - ModelMaskFluxes)**2 / sigma2)

def InterpolateZModel(Z, ModelArray):

    # Define initial Metallicity Array
    ZArray = np.array([Model.Metallicity for Model in ModelArray])
    # Generate Model based on Metallicity
    Model = ModelArray[np.where(ZArray==Z)][0][0] if Z in ZArray else SpectrumModel.InterpolateModels(Z,ModelArray, Linear=True)

    return Model

##### METHODS ##################################################################

def GenerateSteidelMask(Wavelengths):
    ContinuumWindows = Table.read('Parameter Files/Windows/SteidelWindows.txt', format='ascii')
    WindowWlMins = ContinuumWindows['wlmin']
    WindowWlMaxs = ContinuumWindows['wlmax']

    Mask = np.array([np.any((WindowWlMins <= wl) & (wl <= WindowWlMaxs)) for wl in Wavelengths])

    return Mask

##### MAIN METHOD ##############################################################

def FullFitting(stack_path, stack_id, pert_id, model_type, ResultsPath, Verbose=True): #m_id, p_id, model_type, Verbose=True):

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
    #StackSpline = StackModel.DetermineContinuumSpline()

    # Assign Values to easier to reference variables
    Flux = StackModel.Fluxes
    Wavelengths = StackModel.Wavelengths
    Errors = StackModel.Errors

    # Define Steidel Mask
    Mask = GenerateSteidelMask(Wavelengths)
    ModelMask = GenerateSteidelMask(ModelArray[0].Wavelengths)

    # Run Nested Sampling
    with Pool() as pool:
        sampler = dynesty.NestedSampler(LogLikelihood, PriorTransform, ndim=4, logl_args=(Flux, Wavelengths, Errors, ModelArray, Mask, ModelMask), ptform_args=(ZMin, ZMax), pool=pool, queue_size=7)
        sampler.run_nested(print_progress=Verbose)

    # Collect Results + Define Weights
    Results = sampler.results
    Weights = np.exp(Results['logwt'] - Results['logz'][-1])

    # Determine Parameters
    metallicity = dynesty.utils.quantile(x=Results['samples'][:,0], q=[0.16, 0.50, 0.84], weights=Weights)
    Av = dynesty.utils.quantile(x=Results['samples'][:,1], q=[0.16, 0.50, 0.84], weights=Weights)
    B = dynesty.utils.quantile(x=Results['samples'][:,2], q=[0.16, 0.50, 0.84], weights=Weights)
    Delta = dynesty.utils.quantile(x=Results['samples'][:,3], q=[0.16, 0.50, 0.84], weights=Weights)

    # Get Posteriors
    metallicity_samples = dynesty.utils.resample_equal(Results['samples'][:,0], Weights)
    Av_samples = dynesty.utils.resample_equal(Results['samples'][:,1], Weights)
    B_samples = dynesty.utils.resample_equal(Results['samples'][:,2], Weights)
    Delta_samples = dynesty.utils.resample_equal(Results['samples'][:,3], Weights)

    # Output if Verbose
    if Verbose == True: print(metallicity)

    # Record Parameter Data
    with open(f'{ResultsPath}/{ModelPath}/{stack_id}Results.txt', 'a') as file:
        file.write(f'{pert_id}, {metallicity[0]}, {metallicity[1]}, {metallicity[2]}\n')
    with open(f'{ResultsPath}/{ModelPath}/{stack_id}_DustParameters.txt', 'a') as file:
        file.write(f'{pert_id}, Av, {Av[0]}, {Av[1]}, {Av[2]}, B, {B[0]}, {B[1]}, {B[2]}, Delta, {Delta[0]}, {Delta[1]}, {Delta[2]}\n')

    PosteriorTable = Table([metallicity_samples, Av_samples, B_samples, Delta_samples], names=['Z', 'Av', 'B', 'del'])
    PosteriorTable.write(f'{ResultsPath}/SampleDistributions/{stack_id}/p{pert_id}_Posteriors.dat', format='ascii', overwrite=True)
