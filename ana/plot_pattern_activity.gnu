#!/usr/bin/gnuplot

# This is an example to plot the pattern acitivity of the rf2 phase

set border 3
set xtics nomirror 
set ytics nomirror 
set key samplen 1

# Set to your location
datadir = '/home/zenke/data/sim'

set xlabel 'Time (s)'
set ylabel 'Population rate [Hz]'
plot for [i=1:4] sprintf('%s/rf2.0.pact',datadir) using 1:2*i w l t sprintf("Pattern %i",i)
