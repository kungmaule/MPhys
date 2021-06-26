# Relevant Imports
import numpy as np
import sys, os
from astropy.table import Table
from matplotlib import rc
import dynesty as dyn

##### METHODS ##################################################################

def SN(path):
    PathData = Table.read(path, format='ascii')
    Fluxes = PathData['flam']
    Errors = PathData['flam_err']
    return np.mean(Fluxes/Errors)

##### MAIN #####################################################################

if __name__=="__main__":

    if len(sys.argv) == 1:
        print(f'Usage: python {sys.argv[0]} <Root> <Model>')
        quit()
    else:
        Root = sys.argv[1]
        Model = sys.argv[2]

    # Assumes two model types
    ModelType = 'S99' if Model == 'S' else 'BPASS'

    Num = 7 if (Root == 'MassStacks' or Root == 'ContinuumFit_MassStacks') else 14

    if Num == 7:
        Stacks = np.array(['m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7'])
    else:
        Stacks = np.array(['m1sfr1', 'm1sfr2', 'm2sfr1', 'm2sfr2', 'm3sfr1', 'm3sfr2', 'm4sfr1', 'm4sfr2', 'm5sfr1', 'm5sfr2', 'm6sfr1','m6sfr2', 'm7sfr1', 'm7sfr2'])

    # Set up Arrays
    MedianMass = np.empty(Num)
    LowMErr = np.empty(Num)
    HighMErr = np.empty(Num)

    MedianSFR = np.empty(Num)
    LowSFRErr = np.empty(Num)
    HighSFRErr = np.empty(Num)

    MedianSN = np.empty(Num)
    LowSNErr = np.empty(Num)
    HighSNErr = np.empty(Num)
    MinSN = np.empty(Num)

    MedianZ = np.empty(Num)
    LowZErr = np.empty(Num)
    HighZErr = np.empty(Num)

    for count, Stack in enumerate(Stacks):

        if Num == 7:
            DISPATH = f'{Root}/Results/Distributions/{ModelType}/{Stack}_distributions.dat'
        else:
            DISPATH = f'{Root}/Results/{ModelType}Distributions/{Stack}_disributions.dat'

        # Determine Masses
        SFRZM = Table.read(DISPATH, format='ascii')
        Metallicities = SFRZM['lz']
        SFR = SFRZM['lsfr']
        Masses = SFRZM['lmass']

        # Determine Medians and Errors
        TableMetallicities = np.percentile(Metallicities, [16, 50, 84])
        TableMetallicitiesErrors = np.diff(TableMetallicities)
        LowZErr[count], MedianZ[count], HighZErr[count] = TableMetallicitiesErrors[0], TableMetallicities[1], TableMetallicitiesErrors[1]

        TableMasses = np.percentile(Masses, [16, 50, 84])
        TableMassesErrors = np.diff(TableMasses)
        LowMErr[count], MedianMass[count], HighMErr[count] = TableMassesErrors[0], TableMasses[1], TableMassesErrors[1]

        TableSFR = np.percentile(SFR, [16, 50, 84])
        TableSFRErrors = np.diff(TableSFR)
        LowSFRErr[count], MedianSFR[count], HighSFRErr[count] = TableSFRErrors[0], TableSFR[1], TableSFRErrors[1]

        # Determine Average SN
        FileDirectory = f'./{Root}/Stacks/{Stack}'
        FilePaths = [os.path.join(FileDirectory, name) for name in os.listdir(FileDirectory)]

        AverageSN = np.array([SN(path) for path in FilePaths])
        TableSN = np.percentile(SFR, [16, 50, 84])
        TableSNErrors = np.diff(TableSN)
        LowSNErr[count], MedianSN[count], HighSNErr[count] = TableSNErrors[0], np.median(AverageSN), TableSNErrors[1]
        MinSN[count] = np.min(AverageSN)

        print(f'{count+1}/{len(Stacks)} - Done', end='\r')
    print(' -=-=-=-=-=- LOAD IN COMPLETE -=-=-=-=-=-')

    with open(f'TableData/{Root}_{ModelType}.txt', 'a') as file:

        file.write(f'Stack     S/N       min       Mass      Min      Max      SFR      Min      Max       Z          Min        Max\n')
        for i in range(Num):
            ID = Stacks[i]

            M = MedianMass[i]
            LM = LowMErr[i]
            HM = HighMErr[i]

            Z = MedianZ[i]
            LZ = LowZErr[i]
            HZ = HighZErr[i]

            Sf = MedianSFR[i]
            LSf = LowSFRErr[i]
            HSf = HighSFRErr[i]

            Sn = MedianSN[i]
            mSN = MinSN[i]

            if M < 10:
                file.write(f'{ID:6}    {Sn:.3f}    {mSN:.3f}    {M:.4f}    {LM:.3f}    {HM:.3f}    {Sf:.3f}    {LSf:.3f}    {HSf:.3f}    {Z:.5f}    {LZ:.5f}    {HZ:.5f}\n')
            else:
                file.write(f'{ID:6}    {Sn:.3f}    {mSN:.3f}    {M:.3f}    {LM:.3f}    {HM:.3f}    {Sf:.3f}    {LSf:.3f}    {HSf:.3f}    {Z:.5f}    {LZ:.5f}    {HZ:.5f}\n')
