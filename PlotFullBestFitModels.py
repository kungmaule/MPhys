# Relevant Imports
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib import rc
import dynesty as dyn
from SpectrumModel import SpectrumModel, Spline
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)

##### METHODS ##################################################################

def InterpolateZModel(Z, ModelArray):

    # Define initial Metallicity Array
    ZArray = np.array([Model.Metallicity for Model in ModelArray])
    # Generate Model based on Metallicity
    Model = ModelArray[np.where(ZArray==Z)][0][0] if Z in ZArray else SpectrumModel.InterpolateModels(Z,ModelArray, Linear=True)

    return Model

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


def GenerateSteidelMask(Wavelengths):
    ContinuumWindows = Table.read('Parameter Files/Windows/SteidelWindows.txt', format='ascii')
    WindowWlMins = ContinuumWindows['wlmin']
    WindowWlMaxs = ContinuumWindows['wlmax']

    Mask = np.array([np.any((WindowWlMins <= wl) & (wl <= WindowWlMaxs)) for wl in Wavelengths])

    return Mask
##### MAIN #####################################################################

if __name__=="__main__":

    Root = sys.argv[1]
    model = sys.argv[2]

    # Generate Model Array
    if model == 'S':
        ZArray = np.array([0.001, 0.002, 0.008, 0.014, 0.040])
        FWHM, Res = 3.0, 0.4 # 0.4Å Resolution
        folderpath = r"./S99"
        ModelPath = 'S99'
    elif model == 'B':
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

    Limits = [(0.5e-19,1.5e-19), (0.5e-19, 2.0e-19), (0.5e-19, 2.0e-19), (0.5e-19, 2.5e-19), (0.5e-19, 2.5e-19), (0.5e-19, 3.0e-19), (0.5e-19, 3.0e-19)]

    for stack in range(1,8):

        fig, ax = plt.subplots(figsize=(8,6))

        # read in unperturbed stack
        StackModel = SpectrumModel.ReadInTestData(f'{Root}/Stacks/m{stack}/vandels-m{stack}-p0-allz.dat')
        Flux = StackModel.Fluxes
        Wavelengths = StackModel.Wavelengths
        Errors = StackModel.Errors

        ax.plot(Wavelengths, Flux, color='black', lw=.5, label='Stack', zorder=1)
        #ax.plot(Wavelengths, Errors, color='red', lw=.5, label='Errors', zorder=0)

        # Read in Metallicities
        StackData = Table.read(f'{Root}/Results/Distributions/m{stack}_disributions.dat', format='ascii')
        StacklogZDist = StackData['lz']

        # Read in Dust Parameters
        Av = np.empty(0)
        B = np.empty(0)
        d = np.empty(0)
        with open(f'{Root}/Results/{ModelPath}/m{stack}_DustParameters.txt', 'r') as file:
            for count, line in enumerate(file):
                Tokens = line.split(', ')
                Av = np.append(Av, float(Tokens[3]))
                B = np.append(B, float(Tokens[7]))
                d = np.append(d, float(Tokens[11]))

        # Generate Mask
        Mask = ((ModelArray[0].Wavelengths >= Wavelengths[0]) & (ModelArray[0].Wavelengths <= Wavelengths[-1]))
        SteidelDataMask = GenerateSteidelMask(Wavelengths)
        SteidelModelMask = GenerateSteidelMask(ModelArray[0].Wavelengths)

        # Plot unperturbed Stack
        for count, logZ in enumerate(StacklogZDist):
            Z = np.power(10, logZ)
            if Z < 0.001:
                print(Z)
                print(logZ)
            Z_Model = InterpolateZModel(Z, ModelArray)
            A_Lambda = attenuation_salim_2018(wl=Z_Model.Wavelengths/1.e4, av=Av[count], B=B[count], delta=d[count])
            ModelFluxes = Z_Model.Fluxes * np.power(10, -0.4*A_Lambda)
            Beta = np.sum(ModelFluxes[SteidelModelMask]*Flux[SteidelDataMask]/(np.power(Errors[SteidelDataMask],2))) / np.sum(np.power(ModelFluxes[SteidelModelMask]/Errors[SteidelDataMask],2))
            ModelFluxes *= Beta
            # Plot
            if count == 0:
                UnpWl, UnpFlux = Z_Model.Wavelengths[Mask], ModelFluxes[Mask]
                BinVar = np.array(ModelFluxes[Mask])
            else:
                ax.plot(Z_Model.Wavelengths[Mask], ModelFluxes[Mask], color='#3C6997', alpha=0.1, zorder=0)
                BinVar = np.vstack((BinVar, ModelFluxes[Mask]))
        Variance_Bins = np.empty(0)
        for wl in range(BinVar.shape[1]):
            Variance_Bins = np.append(Variance_Bins, np.std(BinVar[:,wl]))
        print(np.mean(Variance_Bins))
        ax.plot(UnpWl, UnpFlux, color='#D62839', label='Unperturbed Best Fit Model', zorder=5, lw=.75)
        ax.set_xlim(1250, 1950)
        ax.set_ylim(Limits[stack-1])
        ax.set_title(f'm{stack}: Plot of All Fits')
        ax.set_ylabel('Flux / $\\textup{erg} \\textup{ s}^{-1} \\textup{cm}^{-2} \\textup{Å}^{-1}$')
        ax.set_xlabel('Wavelength / Å')
        ax.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0))
        plt.tight_layout()
        plt.savefig(f'{Root}/Plots/{ModelPath}/BestFitModels/BestFittingModels_m{stack}.png')
