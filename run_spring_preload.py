import numpy as np
from numpy import random
import time as tme
from multiprocessing import Pool
from simulation3 import simulation
mode = 'springControl'

###############################################################################
# START
#initialization
#
###############################################################################
option = 'const' #choose from const, poly, step(experimental)
name = 'please rename'
path_to_results_directory = 'please rename or leave blank'

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
n_iterations_per_simulation = int(5e4)
beta = 2.
kappa = 15.
stiffness_of_drag_spring = 5.
k_on = 10.
neighborhood_criterion = 0.01
start_time = 0.
distance_between_binding_sites = 2
random.seed(121155)


#scan for drag velocity between
step_min_val = 0
#and
step_max_val = 15
#number of steps in between
step_n_jumps = 5

repetitions_with_same_parameters = 1																																					#|

#configure mulitprocessing TO BE USED WITH CAUTION ----> RAM OVERFLOW
n_cores = 8
use_multiprocessing = False

#//////////////////////////////////////////////////////////////////////////////
# END
#initialization
#//////////////////////////////////////////////////////////////////////////////



#==============================================================================
#SNIPPET
#MULTIPLE RUN, GIVEN DRAG VELOCITY
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
                        th = neighborhood_criterion,
                        d = distance_between_binding_sites,
                        k_pull = stiffness_of_drag_spring)

    velocities = [10,0,-10,0]

    sim.add_run(v_pull=velocities[0])
    s, p, pos, pos_pull = sim.start_run(1)
    sim.add_run(option='loop', v_pull=velocities[1], s=s, p=p, pos_filament=pos, pos_pull=pos_pull)
    s, p, pos, pos_pull = sim.start_run(2, v_pull=velocities[1], s=s, p=p, pos_filament=pos, pos_pull=pos_pull)
    sim.add_run(option='loop', v_pull=velocities[2], s=s, p=p, pos_filament=pos, pos_pull=pos_pull)
    s, p, pos, pos_pull = sim.start_run(3, v_pull=velocities[2], s=s, p=p, pos_filament=pos, pos_pull=pos_pull)
    sim.add_run(option='loop', v_pull=velocities[3], s=s, p=p, pos_filament=pos, pos_pull=pos_pull)
    s, p, pos, pos_pull = sim.start_run(4, v_pull=velocities[3], s=s, p=p, pos_filament=pos, pos_pull=pos_pull)

    n=[1,2,3,4]

    sim.plot_pos(n)
    sim.plot_p(n)
    sim.plot_f(n)
###########################################