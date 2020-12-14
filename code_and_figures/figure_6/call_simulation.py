"""
Iteratively call simulation
"""
import numpy as np
import scipy.io as sio
from Simulation_STP_distparams import runsim
import sys


perc = 8
filename = 'Sim_160000_05/PRR_sim_' + str(perc) + 'perc.mat'

for sim_type in range(4):
    # With free membrane neuron
    spkt, VmP, GI, GE = runsim(filename, sim_type=sim_type)
    file_name = 'Sim_160000_05/New_Neuron_sim_' + str(sim_type) + '_perc_' + str(perc) + '.mat'
    sio.savemat(file_name, {'spkt': spkt, 'VmP': VmP, 'GI': GI, 'GE': GE})

    # # Without free membrane neuron (only spiking)
    # spkt = runsim(filename, sim_type=sim_type)
    # file_name = 'Sim_160000_05/New_Neuron_sim_' + str(sim_type) + '_perc_' + str(perc) + '.mat'
    # sio.savemat(file_name, {'spkt': spkt})
