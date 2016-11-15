#
# $Id: figq.plt 3.38.2.6 1992/11/14 02:25:21 woo Exp $
#
set output 'plot.pdf'
#set terminal postscript eps dashed
set terminal pdfcairo size 6,5 font "Helvetica,12" linewidth 1 rounded fontscale 0.6
analytic = 'vsol1.lst';
numeric = 'shocktube_sol_num.txt';
#
#  Multiplot zur Darstellung der Shocktube-Ergebnisse
#
#
set grid
#set title "Shock-Tube Problem"
#
# Set line style
#
set style line 1 lt rgb "#A00000" lw 1 pt 7
set style line 2 lt rgb "#00A000" lw 1 pt 2
#
set origin 0.0,0.0
set label "shocktube solution at time t = 0.228" at 0.55,2.1
set multiplot
set xrange[.0:1.0]
set key above
#
set size 0.45,0.45
set origin 0.025,0.475
set ylabel "Velocity"
plot analytic u 2:3 title "analytic" ls 1, numeric u 2:3 title "numeric" ls 2
set title
set nolabel
#
set size 0.45,0.45
set origin 0.5,0.475
set nolabel
set ylabel "Density"
plot analytic u 2:4 notitle ls 1, numeric u 2:4 notitle ls 2
#
set size 0.45,0.45
set origin 0.025,0.025
set nolabel
set ylabel "Temperature"
set xlabel "X-Axis"
plot analytic u 2:5 notitle ls 1, numeric u 2:5 notitle ls 2
#
set size 0.45,0.45
set origin 0.5,0.025
set nolabel
set ylabel "Pressure"
set xlabel "X-Axis"
plot analytic u 2:6 notitle ls 1, numeric u 2:6 notitle ls 2
#
set nomultiplot
# undo what we have done above
set title
set autoscale x
set xtics
