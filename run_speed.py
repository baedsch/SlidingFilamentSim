import numpy as np
from numpy import random
import time as tme
from multiprocessing import Pool
from simulation3 import simulation
mode = 'vControl'

###############################################################################
# START
#initialization
#
###############################################################################
option = 'const' #choose from const, poly, step
name = 'please rename'
path_to_results_directory = 'please fill carefully or leave the string blank'

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
neighborhood_criterion = 0.01
start_time = 0.
distance_between_binding_sites = 2.
random.seed(121155)

############################################
#####CASE OPTION CONST######################
#scan for load between
step_min_val = 0
#and
step_max_val = 1.2
#number of steps in between
step_n_jumps = 2

#####CASE OPTION POLY#######################
coeff_of_velocity_polynomial = [[5, 0.], [10, 0.]]
############################################

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
#MULTIPLE RUN, GIVEN VELOCITY
#==============================================================================
if __name__ == '__main__':
    sim = simulation(   mode,
                        n_iterations_per_simulation,
                        n_heads,
                        name = name,
                        path = path_to_results_directory,
                        option = 'const',
                        s_store = s_store,
                        f_store = f_store,
                        f_sum_store = sum_f_store,
                        p_sum_store = sum_p_store,
                        pos_store = pos_store,
                        writeText = True,
                        bta = beta,
                        k = kappa,
                        k_on = k_on,
                        th = neighborhood_criterion,
                        d_t = delta_t,
                        d = distance_between_binding_sites)


    if option == 'const':
        velocities = [i * (step_max_val - step_min_val) / step_n_jumps + step_min_val for i in range(step_n_jumps + 1)]
        for r in range(repetitions_with_same_parameters):
            for v in velocities:
                sim.add_run(v=v, n_sim=r)
        n = [i+1 for i in range(len(velocities) * repetitions_with_same_parameters)]


    if option == 'poly':
        for r in range(repetitions_with_same_parameters):
            for v in coeff_of_velocity_polynomial:
                sim.add_run(v_Coeff=v, n_sim=r)
        n = [i+1 for i in range(len(coeff_of_velocity_polynomial) * repetitions_with_same_parameters)]


    if not use_multiprocessing:
        for ni in n:
            sim.start_run(ni)
    else:
        pool = Pool(processes=n_cores)
        pool.map(sim.start_run, n)
        for nr in n:
            sim.start_run(nr)

    for i in n:
        res = sim.average_norm_force_single(i)
    sim.plot_f_norm__v()
    sim.plot_pos(n)
    sim.plot_p(n)
    sim.plot_f(n)

#==============================================================================

