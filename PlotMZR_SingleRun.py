# Relevant Imports
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib import rc
rc('text', usetex=True)

##### MAIN #####################################################################

if __name__=="__main__":

    Root = sys.argv[1]
    model = sys.argv[2]
    if model.upper() == 'S':
        ModelType = 'S99'
    else:
        ModelType = 'BPASS'

    MassArray = np.empty(shape=(7))
    MetallicityArray = np.empty(shape=(3,7))

    with open(f'{Root}/Data/MassData.txt', 'r') as file:
        for count, line in enumerate(file):
            Tokens = line.split(', ')
            MassArray[count] = float(Tokens[2])

    with open(f'{Root}/Results/{ModelType}/Results.txt', 'r') as file:
        for count, line in enumerate(file):
            Tokens = line.split(', ')
            MetallicityArray[0][count] = np.log10(float(Tokens[2])) - np.log10(float(Tokens[1]))
            MetallicityArray[1][count] = np.log10(float(Tokens[2]))
            MetallicityArray[2][count] = np.log10(float(Tokens[3])) - np.log10(float(Tokens[2]))

    plt.figure()
    plt.errorbar(MassArray, MetallicityArray[1], yerr=(MetallicityArray[0], MetallicityArray[2]), fmt='x', lw=.5, capsize=2)
    a, b = np.polyfit(MassArray, MetallicityArray[1], 1)
    print(a)
    plt.plot(MassArray, a*MassArray+b, lw=.5, linestyle='--')
    plt.ylabel('Log Z')
    plt.xlabel('Log M')
    plt.title('Plot of MZR')
    plt.savefig(f'{Root}/MZR.png')
