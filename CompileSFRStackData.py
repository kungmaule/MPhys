# Imports from libraries
from astropy.table import Table
import numpy as np
import os, sys

##### MAIN #####################################################################

if __name__=="__main__":

    if len(sys.argv) != 2:
        print(f'Usage: python3 {sys.argv[0]} < Mass / SFR >')
        quit()

    Parameter = sys.argv[1]
    Directory = f'./{Parameter}Stacks'
    Splits = ['A', 'B', 'C', 'D', 'E']
    Bins = 7

    # Compile Data
    for stack in range(1, Bins+1):
        for sfr in [1, 2]:
            MassIndex, LowMasses, MedianMasses, HighMasses = [], [], [], []
            SFRIndex, LowSFR, MedianSFR, HighSFR = [], [], [], []
            for count, Split in enumerate(Splits):
                with open(f'{Directory}/MassData/{Split}/sfr{sfr}/m{stack}_MassData.txt', 'r') as file:
                    for line in file:
                        Tokens = line.split(', ')
                        MassIndex.append(int(Tokens[0]))
                        LowMasses.append(float(Tokens[1]))
                        MedianMasses.append(float(Tokens[2]))
                        HighMasses.append(float(Tokens[3]))
                with open(f'{Directory}/SFRData/{Split}/sfr{sfr}/m{stack}_SFRData.txt', 'r') as file:
                    for line in file:
                        Tokens = line.split(', ')
                        SFRIndex.append(int(Tokens[0]))
                        LowSFR.append(float(Tokens[1]))
                        MedianSFR.append(float(Tokens[2]))
                        HighSFR.append(float(Tokens[3]))
            print(f'[ Completion: {stack} / {Bins}]')

            print(len(MassIndex), len(LowMasses), len(MedianMasses), len(HighMasses))
            print(len(SFRIndex), len(LowSFR), len(MedianSFR), len(HighSFR))
            MassTable = Table([MassIndex, LowMasses, MedianMasses, HighMasses], names=['id', 'logMass_Low', 'logMass', 'logMass_High'])
            MassTable.write(f'{Directory}/Results/MassData/m{stack}_sfr{sfr}_MassData.dat', format='ascii', overwrite=True)
            SFRTable = Table([SFRIndex, LowSFR, MedianSFR, HighSFR], names=['id', 'logSFR_low', 'logSFR', 'logSFR_High'])
            SFRTable.write(f'{Directory}/Results/SFRData/m{stack}_sfr{sfr}_SFRData.dat', format='ascii', overwrite=True)
