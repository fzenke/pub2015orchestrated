#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from auryntools import *
from sklearn.decomposition import PCA

# This example performs PCA on the binned spiking activity of excitatory
# cells in the network.
# https://github.com/fzenke/pub2015orchestrated/
# using BinarySpikeMonitors 
# Moreover it assumes that your auryntools is in your PYTHONPATH
# and you have scikitlearn installed

num_mpi_ranks = 4 # the number of sims you used in parallel
datadir = "/home/zenke/data/sim" # Set this to your data path

spkfiles  = ["%s/rf1.%i.e.spk"%(datadir,i) for i in range(num_mpi_ranks)]
sfo = AurynBinarySpikeView(spkfiles)

time_range = 100
bin_size = 100e-3

print("Crunching file ...")
tm = sfo.t_max 
data = sfo.time_binned_spike_counts(tm-time_range,tm,bin_size=bin_size, max_neuron_id=4096)

print("Analyzing ...")
pca = PCA(n_components=5)
pca.fit(data)

print("Plotting ...")
time = np.linspace(tm-time_range,tm,data.shape[0])
plt.plot(time, pca.transform(data))
plt.xlabel("Time (s)")
plt.show()
