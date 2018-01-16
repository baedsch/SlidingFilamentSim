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
name = 'final'
path_to_results_directory = 'please rename or leave blank'

#store data in ram / write them to text files?
s_store = False
p_store = True
f_store = False
sum_f_store = True
sum_p_store = False
pos_store = True
pos_pull_store = True
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
distance_between_binding_sites = 2.
random.seed(121155)

repetitions_with_same_parameters = 5																					#|

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
                        pos_pull_store = pos_pull_store,
                        writeText = writeText,
                        bta = beta,
                        k = kappa,
                        k_on = k_on,
                        th = neighborhood_criterion,
                        d = distance_between_binding_sites,
                        k_pull = stiffness_of_drag_spring)

    for i in range(repetitions_with_same_parameters):
        velocities = [10,0,-10,0]
        t_end  = 0.
        sim.add_run(v_pull=velocities[0], t0=0.)
        s, p, t_end, pos, pos_pull = sim.start_run(1 + 4*i)
        sim.add_run(option='loop', v_pull=velocities[1], s=s, p=p, pos_filament=pos, pos_pull=pos_pull, t0=t_end)
        s, p, t_end, pos, pos_pull = sim.start_run(2 + 4*i, v_pull=velocities[1], s=s, p=p, pos_filament=pos, pos_pull=pos_pull)
        #more steps to ensure that puller passes 0
        n_iterations_per_simulation_extended = int(1.5 * n_iterations_per_simulation)
        sim.add_run(n_steps=n_iterations_per_simulation_extended, option='loopCatchNegative', v_pull=velocities[2], s=s, p=p, pos_filament=pos, pos_pull=pos_pull, t0=t_end)
        s, p, t_end, pos, pos_pull = sim.start_run(3 + 4*i, v_pull=velocities[2], s=s, p=p, pos_filament=pos, pos_pull=pos_pull)
        sim.add_run(option='loop', v_pull=velocities[3], s=s, p=p, pos_filament=pos, pos_pull=pos_pull, t0=t_end)
        s, p, t_end, pos, pos_pull = sim.start_run(4 + 4*i, v_pull=velocities[3], s=s, p=p, pos_filament=pos, pos_pull=pos_pull)

    n=list(range(4 * repetitions_with_same_parameters + 1))[1:]

    sim.plot_pos(n)
    sim.plot_pos_pull(n)
    sim.plot_p(n)
    sim.plot_f(n)
###########################################