# Numerical / System Imports
from astropy.table import Table
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
import numpy as np
from spectres import spectres
from scipy.interpolate import interp1d
import os

# Graphing Imports
from matplotlib import colors, cm, gridspec, rc
import matplotlib.pyplot as plt
#from matplotlib import gridspec
from matplotlib import rc
rc('text', usetex=True)

##### CLASS INSPECIFIC FUNCTIONS ###############################################

def GenerateColours(Seq=False):
    CArr = np.empty(0, dtype='str')
    if Seq == True:
        # Sequential Case
        cmap = cm.get_cmap('plasma', 7)
        for i in range(1, cmap.N-1):
            rgba = cmap(i)
            col = colors.rgb2hex(rgba)
            CArr = np.append(CArr, col)
    else:
        CArr = np.array(['#0A0A0A', '#D62839', '#3CD24C', '#00916E', '#FF9F1C', '#90137E', '#3C6997'])
    return CArr

################################################################################
                              ### CLASSES ###
################################################################################

class SpectrumModel(object):

    def __init__(self, Fluxes, Wavelengths, Errors=None, Metallicity=None):
        self.Metallicity = Metallicity
        self.Fluxes = Fluxes
        self.Wavelengths = Wavelengths
        self.Errors = Errors

##### Self Methods #############################################################

    def ConvolveModel(self, FWHM=3.0, PixelRes=0.4):
        # Convolve the fluxes with a 1DGaussian
        Sigma = (FWHM / np.sqrt(8*np.log(2))) / PixelRes
        Kernel = Gaussian1DKernel(stddev=Sigma)
        self.Fluxes = convolve(np.array(self.Fluxes), Kernel)

    def ConvolveErrors(self, FWHM=3.0, PixelRes=0.4):
        # Convolve the errors with a 1DGaussian
        Sigma = (FWHM / np.sqrt(8*np.log(2))) / PixelRes
        Kernel = Gaussian1DKernel(stddev=Sigma)
        self.Errors = convolve(np.array(self.Errors), Kernel)

    def ResampleRelativeFluxToGrid(self, WavelengthGrid):
        # Rescale fluxes to the new wavelength grid
        self.RelativeFluxes = spectres(WavelengthGrid, self.Wavelengths, self.RelativeFluxes, fill=12345)
        # If Errors, rescale those too
        if self.Errors is not None:
            self.Errors = spectres(WavelengthGrid, self.Wavelengths, self.Errors)
        # Set wavelengths to equal Wavelength grid
        self.Wavelengths = WavelengthGrid

    def ResampleFluxToGrid(self, WavelengthGrid):
        # Rescale fluxes to the new wavelength grid
        self.Fluxes = spectres(WavelengthGrid, self.Wavelengths, self.Fluxes)
        # If Errors, rescale those too
        if self.Errors is not None:
            self.Errors = spectres(WavelengthGrid, self.Wavelengths, self.Errors)
        # Set wavelengths to equal Wavelength grid
        self.Wavelengths = WavelengthGrid

    def DetermineContinuumSpline(self):

        # Determine Window Parameters
        WindowWlMins = np.array([1274.5, 1348.0, 1441.5, 1586.0, 1678.0, 1754.5, 1800.0, 1875.0, 1968.0])
        WindowWlMaxs = np.array([1278.0, 1351.0, 1444.5, 1590.0, 1684.0, 1758.0, 1804.0, 1880.0, 1972.0])
        # Removed :  1485.5, 1488.5,
        # Determine Loop Parameters and Arrays
        WindowNumber = len(WindowWlMins)
        ContinuumFlux = np.empty(WindowNumber)
        WindowWlMeans = np.empty(WindowNumber)

        # Iterate over the windows
        for window in range(WindowNumber):

            RixMask = ((self.Wavelengths >= WindowWlMins[window]) & (self.Wavelengths <= WindowWlMaxs[window]))
            ContinuumFlux[window] = np.median(self.Fluxes[RixMask])
            WindowWlMeans[window] = np.mean([WindowWlMins[window], WindowWlMaxs[window]])

        # Set arrays for the relevant Wavelengths and Fluxes for spline
        SplineMask = ((self.Wavelengths >= WindowWlMeans[0]) & (self.Wavelengths <= WindowWlMeans[-1]))
        SplineWl = self.Wavelengths[SplineMask]
        RelevantFluxes = self.Fluxes[SplineMask]

        # Perform Spline Interpolation
        ContinuumSpline = interp1d(WindowWlMeans, ContinuumFlux, kind='cubic')
        SplineFluxes = ContinuumSpline(SplineWl)

        # Determine Relative Flux Values
        RelativeFluxes = RelevantFluxes / SplineFluxes
        if self.Errors is not None:
            self.Errors = (self.Errors[SplineMask]/SplineFluxes)

        # Return Spline Model
        return Spline(RelevantFluxes, SplineWl, self.Errors, SplineFluxes, RelativeFluxes, Metallicity=self.Metallicity)

##### Static Methods ###########################################################

    def SolarMetallicity(Z):
        return (Z/0.0142)

    def AbsoluteMetallicity(Z):
        return (Z*0.0142)

    def GenerateErrors(Fluxes, SigToNoise=15.0):
        return Fluxes / SigToNoise

    def PerturbFluxes(Fluxes, Errors):
        # Perturb the fluxes within the error values
        PerturbedFluxes = Fluxes + np.random.normal(loc=0, scale=Errors)
        return PerturbedFluxes

    def InterpolateModels(Z, ModelArray, Linear=True):
        # INSERT COMMENT DESCRIBING FUNCTIONS
        InterpolationLimits = np.array([ModelArray[0].Metallicity, ModelArray[-1].Metallicity])

        # Check to see if Z is within range
        if Z > InterpolationLimits[1] or Z < InterpolationLimits[0]:
            print("%f is out of interpolation range." % Z)
            print("Valid metallicity range is %f < Z < %f" % (InterpolationLimits[0], InterpolationLimits[1]))
            quit()

        # Locate the two boundary models.
        for i in range(ModelArray.shape[0]-1):
            if Z > (ModelArray[i].Metallicity) and Z < (ModelArray[i+1].Metallicity):
                LowerModel = ModelArray[i]
                UpperModel = ModelArray[i+1]

        # Set new model arrays
        Wavelengths = LowerModel.Wavelengths
        IntpFluxes = np.empty(Wavelengths.shape[0])

        # Set metallicity limits
        LowerMetallicity = np.log10(LowerModel.Metallicity) # IF NOT ALREADY LOG
        UpperMetallicity = np.log10(UpperModel.Metallicity) # IF NOT ALREADY LOG

        # Iterate over all wavelength values
        ZFrac = (np.log10(Z)-LowerMetallicity)/(UpperMetallicity-LowerMetallicity)
        FluxDifference = UpperModel.Fluxes - LowerModel.Fluxes
        IntpFluxes = LowerModel.Fluxes + (ZFrac * FluxDifference)

        # Return Model with given parameters
        return SpectrumModel(IntpFluxes, Wavelengths, Metallicity=Z)

    ##### DATA READ IN METHODS #

    def ReadInStackData(Path):
        # Read in the data from the stack, and build it into a model
        StackData = Table.read(Path, format='ascii')
        Wavelengths = StackData['Wavelength']
        Fluxes = StackData['Composite Flux']
        Error = StackData['Error Flux']
        return SpectrumModel(Fluxes, Wavelengths, Errors=Error)

    def ReadInTestData(Path):
        StackData = Table.read(Path, format='ascii')
        Wavelengths = StackData['wl']
        Fluxes = StackData['flam']
        Error = StackData['flam_err']
        return SpectrumModel(Fluxes, Wavelengths, Errors=Error)

    def ReadInModelData(Path, Z):
        Data = Table.read(Path, format='ascii')
        #ZDot = (Z/0.0142).round(2)
        Fluxes = Data['f_tot']
        Wavelengths = Data['wl']
        return SpectrumModel(Fluxes, Wavelengths, Metallicity=Z)

################################################################################

################################################################################

class Spline(SpectrumModel):

    def __init__(self, Fluxes, Wavelengths, Errors, SplineFluxes, RelativeFluxes, Metallicity=None):
        self.Fluxes = Fluxes
        self.Wavelengths = Wavelengths
        self.Errors = Errors
        self.Metallicity = Metallicity
        self.SplineFluxes = SplineFluxes
        self.RelativeFluxes = RelativeFluxes

##### METHODS  #################################################################

    def ReScaleErrors(self):
        # Rescales error according to the spline
        self.Errors = (self.Errors / self.SplineFluxes)
