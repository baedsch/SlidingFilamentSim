import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
from numpy import random
#~ import random
import math
import os
import time as tme
# in order to be able to use map with KeyWordArgs
from functools import partial
import multiprocessing as mp


import helper as h
from simulation import simulation


#initialization
########################################################################
mode = 'fControl' #choose from ['vControl', 'fControl', ]
name = 'low_force_check_200_force'
s_store = True
p_store = True
f_store = True
sum_f_store = True
pos_store = True
writeText = True

n_heads = int(5e2)
n_steps = int(7.5e1)
#~ loadF = [.1, 0.5, 1., 5., 10., 50., 100., 500., 1000.]
#~ loadF = [.1, 10., 50., 100.]
#~ loadF = [0. for i in range(2)]
loadF = [200. for i in range(2)]
d_t = 5e-3
bta = 1
k = 30.
th = 0.01
time = 0.
d = 1.
random.seed(121155)

#trying all possibilities for initial p vectors
probabilities_p = [[a / 10., round(1.- a / 10., 2)] for a in range(11)]

n = [i+1 for i in range(len(loadF))]

########################################################################

sim = simulation(mode, n_steps, n_heads, name=name, loadF=loadF, s_store=s_store, f_store=f_store, f_sum_store=sum_f_store, pos_store=pos_store, writeText=writeText, bta=bta, k=k, th=th, d_t=d_t, d=d)
#~ for p in probabilities_p:
	#~ sim.add_run(probabilities_p=p)
for f in loadF:
	sim.add_run(loadF=f)
for rn in n:
	sim.start_run(rn)
	sim.sum_up_P(rn)

q = mp.Queue()
#~ job = mp.Process(target=sim.start_run, args=(1,))
#~ job.start()
#~ job.join()
#~ jobs = [mp.Process(target=sim.start_run, args=(i,q)) for i in n]
#~ for j in jobs:
	#~ j.start()
#~ for j in jobs:
	#~ j.join()
#~ sim.create_threadsMP(n)
#~ sim.start_threadsMP()
#~ print sim.sum_P[1]
#~ print sim.t[1]
sim.plot_pos(1)
sim.plot_p(1)
sim.plot_f(1)

print probabilities_p

########### ONE SINGLE RUN SNIPPET
#~ sim.start_run(0)
#~ sim.sum_up_P(0)
#~ sim.plot_f(0)
#~ sim.plot_p(0)
#~ sim.plot_pos(0)
##################################
