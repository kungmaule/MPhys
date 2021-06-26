from astropy.io import fits
import numpy as np
from astropy.cosmology import WMAP9 as cosmo

class VANDELS_Galaxy(object):

    """
    Class to describe a a VANDELS Galaxy and it's data

    Properties:
    *List Properties when added

    Methods:
    *List methods when written*
    """

    def __init__(self, hdu, data, hdr):
        """
        Initialise instance of galaxy with relevant quantities
        """
        #self.name = name
        self.hdu = hdu
        self.data = data
        self.hdr = hdr

        self.Flux = hdu['EXR1D'].data
        self.Flux_Error = hdu['NOISE'].data
        self.Rescaled_Flux = None
        self.Redshift = float(hdr['HIERARCH PND Z'])


        self.wl0 = self.hdr['CRVAL1']
        self.dwl = self.hdr['CDELT1']
        self.naxis = self.Flux.shape[0]
        self.wlf = self.wl0 + (self.dwl*self.naxis)

        # Determine Range
        self.Lamda_Obs = np.arange(self.wl0, self.wlf, self.dwl)
        # Shift to rest frame
        self.Lamda_RestFrame = self.Lamda_Obs / (1.0 + self.Redshift)



    def __str__(self):
        """
        Return Galaxy details in string format
        """

        return ("Observed Wavelength: " + str(Lamda_Obs) + "\nRedshift: " + str(Redshift) + "\nRest Frame Wavelength: " + str(Lamda_RestFrame))

    ########################## METHODS ##########################

    def Convert_Flux_To_Median_Z(self, Z_Median):
        # Determine Luminosity Distance at actual redshift
        L_Distance = cosmo.luminosity_distance(self.Redshift)
        # Determine Luminosity Distance at Median redshift
        L_Median_Redshift = cosmo.luminosity_distance(Z_Median)
        # Determine New Flux
        self.Rescaled_Flux = self.Flux * (L_Distance/L_Median_Redshift)**0.5

    def Add_Parameters(self, Mass, Mass_Error, SFR, SFR_Error):
        self.logMass = Mass
        self.logMass_Error = Mass_Error
        self.logSFR = SFR
        self.logSFR_Error = SFR_Error

    def Set_Name(self, name):
        self.name = name

    ####################### STATIC METHODS ######################

    @staticmethod
    def Generate_Galaxy(path):
        hdu = fits.open(path)
        data = hdu[1].data
        hdr = hdu[0].header
        return VANDELS_Galaxy(hdu, data, hdr)
