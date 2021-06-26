Game of Life : Read Me Document

Please find in the following a description and explanation of my submission for this checkpoint.

USAGE: SIRS.py <Lattice Size> <Run Mode [S, P, C, I]> <Animate(Y/M)>

Some user error handling has been implemented to ensure correct inputs. The inputs are as follows:
- There are various runmodes. 'S' performs a single run for user specified probabilities. 'P' fixes p2=0.5 and generates the p1p3 plane data
  for a user spectified pInterval. 'C' does the same as 'P' but fixes p3=0.5 and iterates over a smaller range of p1. 'I' iterates over an
  immune fraction of user specified interval.
- Picking a lattice size of N will generate an NxN lattice of cells. These will initially be randomly assigned as either immune (2, only if
  chosen to be included), recovered (1), susceptible (0) or infected (-1).
- Inputting "Y" will cause the lattice to be animated and output to the user. Inputting "N" will skip this step.

Some values have been hard-coded.
- Equilibration steps = 100
- StopIfDead will kill the infectivity run if the infection dies. It is only false for S, for the purposes of genreating a plot for an
  absorbing state.
- p2, p3 fixed at 0.5 when needed for 'P' or 'C'. For 'C', p1 iterates between 0.2 and 0.5.
Whilst these values are hard-coded, they are simple to change, allowing the code to be very flexible at modelling various dynamics.

The single run method performs a run at a specific set of p1, p2, p3 specified by a user. Example parameters to run for some states are:
- Dynamic Equilibrium: p1 = p2 = p3 = 0.5
- Absorbing: p1, p2, p3 = 0.2, 0.5, 0.2
- Waves: p1, p2, p3 = 0.5, 0.1, 0.01
This method does not stop if dead so that it can record the data for the behaviour plots.

The Plane method fixes p2=0.5 and iterates over a user-specified interval in p1 and p3, recording the infectivity fraction and variance.
The run for a given parameter set will stop if the infection dies, setting the infectivity and variance to zero. The data is recorded to Data,
from which the plane plot can be produced for both fractional infectivity and variance. Two example plots are given.

The cut mode fixes p2, p3 = 0.5 and iterates over a shorter range of p1. The variance is recorded, and error bars calculated using a resampling
method. The data is recorded into Data. Two example plots are given. Note that error bars may be underestimated by the bootstrapping in
comparison to if multiple simulations were run.

The immunity mode introduces a fourth state into the lattice (2) which is a perma-immunity. These states do not interact with neighbouring states.
They are randomly assigned using np.random.choice at the start. Note that in theory, immunity could also be added to any of the above modes. This mode iterates
over the immune fraction, increasing it each run, recording infectivity. For each run, the simulation actually runs 5 times and the standard error on the
mean is taken to be the error. Data is saved to Data.

All data corresponding to plots can be found in Data/Backup.

Final Results:

- Absorbing State Parameters: p1, p2, p3 = 0.2, 0.5, 0.2 (See Plane for any 0 infectivity spots)
- Dynamic Equilibrium Parameters: p1 = p2 = p3 (e.g. 0.5)
- Waves Parameters: p1, p2, p3 = 0.5, 0.1, 0.01 (see WaveLocations)

- Waves locations from the p1p3 plane: See WaveLocations
- Waves parameters from the location in the cut: ~ p1 = 0.33

- Minimum immune fraction for absorbing state [(0.5, 0.5, 0.5) size 50 lattice]: 0.26
- Minimum immune fraction for absorbing state [(0.8, 0.1, 0.02) size 100 lattice]: 0.35
