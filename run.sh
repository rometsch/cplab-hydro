#!/bin/bash

python3 sim_advection.py

python3 loop_vs_array.py

python3 sim_shocktube.py

gnuplot plot_shocktube.plt
