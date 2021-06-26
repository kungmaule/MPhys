################################################################################

"""             Thomas Stanton - Fit Running Code                            """

################################################################################

# Imports from libraries
import numpy as np
import os, sys

# Import Method
from FullFitting import FullFitting

##### MAIN #####################################################################

if __name__=="__main__":

    # Input Handling
    if len(sys.argv) == 1:
        print(f'Correct Usage: python {sys.argv[0]} <Directory> <Mass Stack ID> (<SFR ID>) <Model ID (S/B)>')
        quit()

    # Assign Variables
    Directory = sys.argv[1]
    xStackID =sys.argv[2]
    Model_ID = sys.argv[3].upper()

    for StackID in ['m1sfr1', 'm1sfr2', 'm2sfr1', 'm2sfr2', 'm3sfr1', 'm3sfr2', 'm4sfr1', 'm4sfr2', 'm5sfr1', 'm5sfr2', 'm6sfr1', 'm6sfr2', 'm7sfr1', 'm7sfr2']:
        # Set up file pathing
        FileDirectory = f'./{Directory}/Stacks/{StackID}'
        ResultsPaths = f'./{Directory}/Results'
        FilePaths = [os.path.join(FileDirectory, name) for name in os.listdir(FileDirectory)]

        # Iterate over all Stacks
        print(f'Fitting {len(FilePaths)} perturbations of the {StackID} stack.')
        for count, stack in enumerate(FilePaths):
            # Fit Metallicity
            FullFitting(stack_path=stack, stack_id=StackID, pert_id=count, model_type=Model_ID, ResultsPath=ResultsPaths, Verbose=False)
            # Update Output
            print(f'[ Stack {count+1}/{len(FilePaths)} fitted. Overall Completion: {100*(count+1)/len(FilePaths)} % ]', end='\r')
        print('\n')
