# Imports from libraries
from astropy.table import Table
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from astropy.stats import sigma_clip
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
import os
import time
import sys, copy

# Imports from Classes
from VANDELS_Gal import VANDELS_Galaxy

##### METHODS ##################################################################

def BinByMass(GalArr, Number):

    SortedGals = sorted(GalArr, key=lambda Gal: Gal.logMass)
    MassStacks = np.split(SortedGals, [151, 303, 455, 606, 758, 910])

    return MassStacks

def Median_Redshift(VANDELS_Galaxies):
    redshifts = np.array([Gal.Redshift for Gal in VANDELS_Galaxies])
    #Z_Array = np.empty(VANDEL_Galaxies.shape[0])
    #for i in range(VANDEL_Galaxies.shape[0]):
        #Z_Array[i] = VANDEL_Galaxies[i].Redshift
    return np.median(redshifts)

def StackSpectra(Stack, m, p):
    # Define Loop parameters
    Initial_Lamda = 1100
    Final_Lamda = 2000
    Resolution = 1
    Wavelengths = np.arange(Initial_Lamda, Final_Lamda, Resolution)

    ConcWavelengths, ConcFluxes = np.empty(0), np.empty(0)
    for Gal in Stack:
        ConcWavelengths = np.concatenate((ConcWavelengths, Gal.Lamda_RestFrame))
        ConcFluxes = np.concatenate((ConcFluxes, Gal.Flux))

    Mask = np.array([np.any((Initial_Lamda <= wl) & (wl <= Final_Lamda)) for wl in ConcWavelengths])
    MaskWavelengths = ConcWavelengths[Mask]
    MaskFluxes = ConcFluxes[Mask]

    Composite_Flux = np.empty(Wavelengths.shape[0])
    Error_Flux = np.empty(Wavelengths.shape[0])

    for count, Pixel in enumerate(Wavelengths):
        Pixel_Flux = ConcFluxes[np.where(np.logical_and(Pixel <= ConcWavelengths, ConcWavelengths <= (Pixel +1)))]
        ClippedPixelFlux = sigma_clip(Pixel_Flux, sigma=3, masked=False)
        Error_Flux[count] = np.std(bootstrap(ClippedPixelFlux, 1000, bootfunc=np.median))
        Composite_Flux[count] = np.median(ClippedPixelFlux)

    StackFile = Table([Wavelengths, Composite_Flux, Error_Flux], names=['wl', 'flam', 'flam_err'])
    StackFile.write(f'../FirstMZRTest/Spectra/m{m}/vandels-m{m}-p{p}-allz.dat', format='ascii', overwrite=True)

##### MAIN METHOD ##############################################################

if __name__=="__main__":

    Set = int(sys.argv[1])

    # Read in Parameter File
    #ParamFile = Table.read('VANDELS_DATA/vandels_tstanton_dr3.dat', format='ascii')
    #SpecNames = ParamFile['spectrum']
    #Redshifts = ParamFile['redshift']
    #Masses = ParamFile['lmass']
    #SFRs = ParamFile['lsfr']

    # Read in Cullen parameter file
    CParamFile = Table.read('params.sedparams', format='ascii')
    SpecNames = CParamFile['id']
    Masses = CParamFile['lmass_fpp']
    SFRs = CParamFile['lsfr_fpp']


    # Hard Coded Errors in Mass, SFR
    M_Err = 0.15
    SFR_Err = 0.30

    # Hard Code Number of Stacks + Perturbations
    StackNumber = 7
    Perturbations = 500

    # Determine Folder path
    VANDELS_Galaxies = np.empty(0)
    print('Reading in Galaxies...')
    for count, path in enumerate(SpecNames):
        Galaxy = VANDELS_Galaxy.Generate_Galaxy(f'VANDELS_DATA/{path}')
        Galaxy.Set_Name(path)
        Galaxy.Add_Parameters(Masses[count], M_Err, SFRs[count], SFR_Err)
        if np.isnan(Galaxy.logMass) or np.isnan(Galaxy.logSFR):
            continue
        VANDELS_Galaxies = np.append(VANDELS_Galaxies, Galaxy)
        print('[Completion: ' + str(round(100*count/len(SpecNames),3)) + ' %]', end='\r')
    print('Galaxy load in completed.')

    # Generate Arrays for storing Values
    MassArray = np.zeros(shape=(StackNumber, Perturbations, 3))

    # Iterate over Perturbations
    if Set == 1:
        pRange = range(0,50,1)
        X = 'A'
    elif Set == 2:
        pRange = range(100, 200, 1)
        X = 'B'
    elif Set == 3:
        pRange = range(200, 300, 1)
        X = 'C'
    elif Set == 4:
        pRange = range(300,400,1)
        X = 'D'
    elif Set == 5:
        pRange = range(400, 500,1)
        X = 'E'

    print(f'Beggining Stacking of set: {X} over {pRange}')
    for pert in pRange:


        if pert == 0:
            # Don't perturb
            StackedGalObjs = BinByMass(VANDELS_Galaxies, StackNumber)
            #with open(f'../FirstMZRTest/Numbers.txt', 'a') as file:
                #for count, Stack in enumerate(StackedGalObjs):
                    #file.write(f'{count+1}: {Stack.shape[0]} Galaxies\n')
        else:
            Perturbed_VANDELS_Galaxies = copy.deepcopy(VANDELS_Galaxies)
            for Gal in Perturbed_VANDELS_Galaxies:
                Gal.logMass = np.random.normal(loc=Gal.logMass, scale=Gal.logMass_Error)

            StackedGalObjs = BinByMass(Perturbed_VANDELS_Galaxies, StackNumber)

        for count, Stack in enumerate(StackedGalObjs):

            # Record Mass Data
            StackMass = np.array([Gal.logMass for Gal in Stack])
            MassArray[count, pert, 0] = np.min(StackMass)
            MassArray[count, pert, 1] = np.median(StackMass)
            MassArray[count, pert, 2] = np.max(StackMass)

            with open(f'../FirstMZRTest/Data/m{count}_MassData.txt', 'a') as file:
                file.write(f'{pert}, {np.min(StackMass)}, {np.median(StackMass)}, {np.max(StackMass)} \n')

            # Determine Median Redshift of the Stack + Normalise Flux
            zMedian = Median_Redshift(Stack)
            for Gal in Stack:
                Gal.Convert_Flux_To_Median_Z(zMedian)

            StackSpectra(Stack, count+1, pert)
            print('[Completion:' + str(round(100*(count+1)/len(StackedGalObjs),2)) + '%]', end='\r')

        print(f' ------ Pertubation {pert} Completed ------\n')
        print('[Completion:' + str(round(100*pert/Perturbations,2)) + '%]', end='\r')
