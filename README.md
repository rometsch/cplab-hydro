# cplab-hydro
hydrodynamics project for my computational physics lab

Run `run.sh` to produce all results of this lab project.
Modify this script to use the right commands for python and gnuplot.

## shocktube problem
A numerical solution of the shocktube problem is performed in `sim_shocktube.py`.
To run the simulation and produce a plot comparing the obtained solution to an
analytical one run
`run.sh`
on a unix system.
You'll need `gnuplot` with the `pdfcairo terminal` and Python3.

## advection
Execute `sim_advection.py` to run differnt methods for linear advection and to get a plot comparing the differnt methods.

## loop vs array
`loop_vs_array.py` is a script to compare execution times of the different methods.
Especially the two different implementations of the 2nd order upwind method are of interest.
The implementation using Numpy array multiplication is about 30-40 times faster than a naive implementation using loops.
