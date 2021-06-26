# Numerical / System Imports
from astropy.table import Table
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
import numpy as np
from spectres import spectres
from scipy.interpolate import interp1d
import os
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy import constants as c
# Graphing Imports
from matplotlib import colors, cm, gridspec, rc
import matplotlib.pyplot as plt
#from matplotlib import gridspec
from matplotlib import rc
rc('text', usetex=True)
# Imports from Classes

from SpectrumModel import SpectrumModel, Spline
from VANDELS_Gal import VANDELS_Galaxy

##### METHODS ##################################################################
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

def C19Metallicity(logMass):
    # Assuming no scatter
    m10 = logMass - 10
    logZperSolar = (0.27*m10) -0.68
    logZ = logZperSolar + np.log10(0.0142)
    return np.power(10, logZ)

def InterpolateZModel(Z, ModelArray):

    # Define initial Metallicity Array
    ZArray = np.array([Model.Metallicity for Model in ModelArray])
    # Generate Model based on Metallicity
    if Z < ModelArray[0].Metallicity:
        Model = ModelArray[0]
    elif Z > ModelArray[-1].Metallicity:
        Model = ModelArray[-1]
    else:
        Model = ModelArray[np.where(ZArray==Z)][0][0] if Z in ZArray else SpectrumModel.InterpolateModels(Z,ModelArray, Linear=True)

    return Model

##### MAIN #####################################################################

if __name__=="__main__":

    # Read in Models
    ZArray = np.array([0.001, 0.002, 0.008, 0.014, 0.040])
    FWHM, Res = 3.0, 0.4
    folderpath = r"./S99"
    filepaths = [os.path.join(folderpath, name) for name in os.listdir(folderpath)]
    ModelArray = np.empty(0)
    ModelWavelengths = np.arange(1000, 2000, 1)
    for count, path in enumerate(filepaths):
        # Read in Data and generate SpectrumModel
        Model = SpectrumModel.ReadInModelData(path, ZArray[count])
        ModelArray = np.append(ModelArray, Model)

    # Read in Data Parameters
    CParamFile = Table.read('StackingMethods/params.sedparams', format='ascii')
    SpecNames = CParamFile['id']
    Masses = CParamFile['lmass_fpp']
    SFRs = CParamFile['lsfr_fpp']
    # Hard Coded Errors in Mass, SFR
    M_Err = 0.15
    SFR_Err = 0.30

    #VANDELS_Galaxies = np.empty(0)
    SimGal_ID = np.empty(0)
    SimGal_Masses = np.empty(0)
    SimGal_Metallicities = np.empty(0)
    SimGal_SFR = np.empty(0)
    SimGal_SN = np.empty(0)
    SimGal_z = np.empty(0)

    # set the cosmology assuming H0=70 km/s/Mpc and Omega_m=0.3
    # this is for calculating luminosity distances
    cosmo = FlatLambdaCDM(70, 0.3)

    print('Reading in Galaxies...')
    GalCount = 0
    for count, path in enumerate(SpecNames):

        # Skip nans
        if np.isnan(Masses[count]) or np.isnan(SFRs[count]):
            continue

        # Read in Individual
        Galaxy = VANDELS_Galaxy.Generate_Galaxy(f'StackingMethods/VANDELS_DATA/{path}')
        Galaxy.Set_Name(path)
        Galaxy.Add_Parameters(Masses[count], M_Err, SFRs[count], SFR_Err)
        GalCount += 1

        # Set up variables
        logMass = Galaxy.logMass
        Z = C19Metallicity(logMass)
        logSFR = Galaxy.logSFR
        ID = GalCount
        Redshift = Galaxy.Redshift

        # Generate S99 Model of correct Metallicity
        InterpolatedModel = InterpolateZModel(Z, ModelArray)

        # Convolve to correct resolution
        InterpolatedModel.ConvolveModel(FWHM=FWHM, PixelRes=Res)
        Lambda_Obvs = Galaxy.Lamda_RestFrame * (1 + Redshift)

        ##################################################################

        """ Scale luminosity into the correct units """

        # Convert Luminosity to relevant units
        wl = InterpolatedModel.Wavelengths
        lum_sol = InterpolatedModel.Fluxes

        # Convert luminsoty to erg/s
        lum = lum_sol * c.L_sun.to(u.erg/u.s).value

        # Determine lum_scale_factor
        lum_scale_factor = np.power(10, logMass-8.0)

        # Scale Luminosity
        lum *= lum_scale_factor

        # Get luminosity distance
        d_l = cosmo.luminosity_distance(z=Redshift).to(u.cm).value

        # Convert to flux density at the redshift
        flux = (lum / (4. * np.pi * d_l**2)) * (1 / (1. + Redshift))

        ##################################################################

        """ Handle Dust Attenuation """

        x = logMass - 10
        av = (2.293 + 1.160*x + 0.256*x*x + 0.209*x*x*x) / 2.5

        # Scatter the value by 0.3 dex
        av += np.random.normal(0.0, scale=0.3)
        if av < 0.:
            av = 0.

        # Convert to microns
        wl_microns = wl / 1.e4
        alam = attenuation_salim_2018(wl=wl_microns, av=av, B=0.0, delta=0.0)

        # Attenuate the Flux
        flux *= np.power(10, -0.3 * alam)

        ###################################################################

        ScaledError = 2.5*Galaxy.Flux_Error
        PerturbedFluxes = spectres(Lambda_Obvs, wl*(1+Redshift), flux, verbose=False) + np.random.normal(loc=0, scale=ScaledError)

        SimulatedModel = SpectrumModel(PerturbedFluxes, Lambda_Obvs, ScaledError, Z)

        SN = np.mean(SimulatedModel.Fluxes/SimulatedModel.Errors)

        SimGal_ID = np.append(SimGal_ID, ID)
        SimGal_Masses = np.append(SimGal_Masses, logMass)
        SimGal_Metallicities = np.append(SimGal_Metallicities, Z)
        SimGal_SFR = np.append(SimGal_SFR, logSFR)
        SimGal_SN = np.append(SimGal_SN, SN)
        SimGal_z = np.append(SimGal_z, Redshift)

        # Save Galaxy Data

        OutputTable = Table([SimulatedModel.Wavelengths, SimulatedModel.Fluxes, SimulatedModel.Errors], names=['wl', 'flam', 'flam_err'])
        OutputTable.write(f'./VANDELS_SIM/Galaxies/Gal_ID{ID}.dat', format='ascii', overwrite=True)

        print('[Completion: ' + str(round(100*count/len(SpecNames),3)) + ' %]', end='\r')
    print('Galaxy load in completed.')

    # Save overall Data
    ParamTable = Table([SimGal_ID, SimGal_Masses, SimGal_SFR, SimGal_Metallicities, SimGal_z], names=['id', 'lmass', 'lsfr', 'z_stellar', 'z'])
    ParamTable.write('./VANDELS_SIM/ParameterTable.dat', format='ascii', overwrite=True)
