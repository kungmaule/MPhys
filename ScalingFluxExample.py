import numpy as np

from astropy.table import Table

from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy import constants as c

import matplotlib.pyplot as plt

# set the cosmology assuming H0=70 km/s/Mpc and Omega_m=0.3
# this is for calculating luminosity distances
cosmo = FlatLambdaCDM(70, 0.3)

# assume the galaxy mass is 10^(10.4) M_sun and a redshift 3.5
lmass_galaxy = 10.4
redshift = 3.5


# loading a model (filepath specific to my laptop):
# in reality you will have generated a models at the correct metallicity value
sb99 = Table.read('/Users/fcullen/work/vandels/sps_models/starburst99/S99-v00-Z001-IMF2.3_100myr.spectrum',
                 format='ascii.commented_header')

# wavelength
wl = sb99['wl']
lum_sol = sb99['f_stellar'] # the defaul units are solar luminosities

# convert to luminosity in ergs/s
lum = lum_sol * c.L_sun.to(u.erg/u.s).value

# the mass of the model is 10^8 M_sun so need to scale up/down appropriately
# this equation holds for all values of lmass_galaxy
lum_scale_factor = np.power(10, lmass_galaxy - 8.0)

# scale the luminosity:
lum *= lum_scale_factor

# now get the luminsoity distance in cm based on the redshift:
d_l = cosmo.luminosity_distance(z=redshift).to(u.cm).value

# now can convert to flux density units erg/s/cm^2/A at the redshift:
flux = (lum / (4. * np.pi * d_l**2)) * (1 / (1. + redshift))

# I forgot something else: dust attenuation!
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

# get the value of attenuation (Av) using the realtion
# from McLure et al. 2018 (https://ui.adsabs.harvard.edu/abs/2018MNRAS.476.3991M/abstract)
# which gives Av for a given mass
x = lmass_galaxy - 10
av = (2.293 + 1.160*x + 0.256*x*x + 0.209*x*x*x) / 2.5

# scatter this value by 0.3 dex:
av += np.random.normal(0.0, scale=0.3)

# but it can't be negative:
if av < 0.:
    av = 0.

# get A(lambda) using the Salim relation assuming the bump size = 0 and the delta parameter = 0
# remember wavelegnth needs to be in microns for this equation
wl_microns = wl / 1.e4
alam = attenuation_salim_2018(wl=wl_microns, av=av, B=0.0, delta=0.0)

# now attenuate the flux
flux *= np.power(10, -0.4 * alam)

# plot and check the flux level is roughly correct, should be around 10^-19/10^-18 erg/s/cm^2/A for
# the vandels galaxies
plt.plot(wl, flux)
plt.show()

# next tasks:
# - convolve to VANDESL resolution
# - shift to observed frame (wl * (1 + redshift))
# - perturb using the actual VANDELS error spectrum
