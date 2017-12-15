import numpy as np
from numpy import random
import time as tme
from multiprocessing import Pool
from simulation3 import simulation
mode = 'fControl'

###############################################################################
# START
#initialization
#
###############################################################################
option = ''
name = 'loadrange__0_90__direct_detachment'
path_to_results_directory = 'please fill carefully or comment leave the string blank'

#store data in ram / write them to text files?
s_store = False
p_store = True
f_store = False
sum_f_store = True
sum_p_store = False
pos_store = True
writeText = True

#most important parameters
n_heads = int(1e2)
n_iterations_per_simulation = int(1e4)
delta_t = 5e-3
beta = 2.
kappa = 10.
k_on = 10.
neighbourhood_criterion = 0.01
start_time = 0.
distance_between_binding_sites = 2.
random.seed(121155)

#scan for load between
step_min_val = 0
#and
step_max_val = 20
#nubber of steps inbetween
step_n_jumps = 4

repetitions_with_same_parameters = 10																																						#|

#configure mulitprocessing TO BE USED WITH CAUTION ----> RAM OVERFLOW
n_cores = 8
use_multiprocessing = False

#//////////////////////////////////////////////////////////////////////////////
# END
#initialization
#
#//////////////////////////////////////////////////////////////////////////////

#==============================================================================
#SNIPPET
#MULTIPLE RUN, GIVEN LOAD
#==============================================================================
if __name__ == '__main__':
    sim = simulation(   mode,
                        n_iterations_per_simulation,
                        n_heads,
                        name = name,
                        path = path_to_results_directory,
                        option = option,
                        s_store = s_store,
                        f_store = f_store,
                        f_sum_store = sum_f_store,
                        pos_store = pos_store,
                        writeText = writeText,
                        bta = beta,
                        k = kappa,
                        k_on = k_on,
                        th = neighbourhood_criterion,
                        d_t = delta_t,
                        d = distance_between_binding_sites)

    loads = [i * (step_max_val - step_min_val) / step_n_jumps + step_min_val for i in range(step_n_jumps + 1)]
    for r in range(repetitions_with_same_parameters):
        for f in loads:
            sim.add_run(loadF=f, n_sim=r)
    n = [i+1 for i in range(len(loads) * repetitions_with_same_parameters)]
    for rn in n:
    	sim.start_run(rn)
    sim.plot_v__f_norm()
    sim.plot_pos(n)
    sim.plot_p(n)
    sim.plot_f(n)
#==============================================================================