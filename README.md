# cplab-hydro
hydrodynamics project for my computational physics lab

## sim
Execute `sim.py` to run differnt methods for linear advection and to get a plot comparing the differnt methods.

## loop vs array
`loop_vs_array.py` is a script to compare execution times of the different methods.
Especially the two different implementations of the 2nd order upwind method are of interest.
The implementation using Numpy array multiplication is about 30-40 times faster than a naive implementation using loops.
