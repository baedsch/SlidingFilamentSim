import numpy as np
from numpy import random
import time as tme
import multiprocessing as mp
from multiprocessing import Pool

from simulation3 import simulation#, Consumer#, varstep_sim

#initialization

#################################################################################################|

mode = 'vControl' #choose from ['vControl', 'fControl', ]
option = 'const' #for fControl choose from xy			
				#for vControl choose from 													
				#						-> poly: specify coefficients						
				#						-> step: specify n_elem, n_jumps, min_val, max_val		
name = 'spring_crosscheck_v=5'													
																								
#store data in ram / write them to text files?													
s_store = False																					
p_store = True																					
f_store = True																					
sum_f_store = True																				
sum_p_store = True																				
pos_store = True																			     	
writeText = False																			
																								
#most important parameters																		
n_heads = int(1e2)																			
n_steps = int(1e3)																				
d_t = 5e-3																						
bta = 2.																						
k = 10.
k_on = 10.																						
th = 0.01																						
t0 = 0.																						
d = 2.
repetitions = 10																		
random.seed(121155)																				
																								
#parameters for fControl																		    
loadF = [10. for i in range(3)]																
																								
#parameters for vControl																		    
	#-> poly option																				
v_Coeff = [[30, 0.] for i in range(3)] #example for constant velocity of 1.										
																								
	#-> step option																				
#step_n_jumps = int(n_steps / 1000)															    
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
#step_n_jumps = [25.,25]																		    
#step_min_val = [0., 0.05]
#step_max_val = [0.05, 0.2]
#step_n_jumps = [25., 25.]																		    
#step_min_val = [0., 0.05]																		
#step_max_val = [0.05, 0.2]

    #-> const option																			
step_n_jumps = 25															    
step_min_val = 0
step_max_val = 1000

  
#step_factors = [1, 3]																		    
if option == 'poly': 
    if not len(step_min_val) == len(step_max_val) and mode == 'vControl': 		                    
        raise ValueError("step_min_val, step_max_val and step_factors must have the same length")																							#|
                                    																
#configure mulitprocessing																		
n_cores = 1 																			
																								
																								
#################################################################################################|

########### MULTIPLE RUN, GIVEN VELOCITY SNIPPET

#sim = simulation(mode,
#					n_steps,
#					n_heads,
#					name = name,
#					option = option,
#					s_store = s_store,
#					f_store = f_store,
#					f_sum_store = sum_f_store,
#					pos_store = pos_store,
#					writeText = writeText,
#					bta = bta,
#					k = k,
#                k_on = k_on,
#					th = th,
#					d_t = d_t,
#					d = d,
#                repetitions = repetitions,
#                step_n_jumps = step_n_jumps,
#                step_min_val = step_min_val,
#                step_max_val = step_max_val)
#velocities = [i * (step_max_val - step_min_val) / step_n_jumps + step_min_val for i in range(step_n_jumps + 1)]
#for v in velocities:
#	sim.add_run(v=v)
#n = [i+1 for i in range(len(velocities))]
#for rn in n:
#	sim.start_run(rn)
#sim.plot_pos(n)
#sim.plot_p(n)
#sim.plot_f(n)
###########################################

########### MULTIPLE RUN, GIVEN VELOCITY RANGE , #STEPS SNIPPET, MP
if __name__ == '__main__':
    sim = simulation(mode,
    					n_steps,
    					n_heads,
    					name = name,
    					option = 'const',
    					s_store = s_store,
    					f_store = f_store,
    					f_sum_store = sum_f_store,
    					p_sum_store = sum_p_store,
    					pos_store = pos_store,
    					writeText = True,
    					bta = bta,
    					k = k,
                    k_on = k_on,
    					th = th,
    					d_t = d_t,
    					d = d,
                    v = 0,
                    repetitions = 1,
                    step_n_jumps = step_n_jumps,
                    step_min_val = step_min_val,
                    step_max_val = step_max_val)
    
    velocities = [i * (step_max_val - step_min_val) / step_n_jumps + step_min_val for i in range(step_n_jumps + 1)]
    
    
    for r in range(repetitions):
        for v in velocities:
            sim.add_run(v=v, n_sim=r)
    n = [i+1 for i in range(len(velocities))]
    pool = Pool(processes=n_cores)
    pool.map(sim.start_run, n)
#    for nr in n: sim.start_run(nr)
#    for i in n:
#        res = sim.average_norm_force_single(i)
    sim.plot_f_norm__v()
#    Processes = []
#    
#    
#    if len(velocities) % n_cores == 0:
#        runs_per_core = len(velocities) / n_cores
#        for i in range(n_cores):
#            Processes.append(mp.Process(target = sim.start_run, args = (n[:runs_per_core],)))
#            n = n[runs_per_core:]
#    else:
#        runs_per_core = len(velocities) // (n_cores - 1) #floor divided
#        for i in range(n_cores - 1):
#            Processes.append(mp.Process(target = sim.start_run, args=(n[:runs_per_core],)))
#            n = n[runs_per_core:]
#        Processes.append(mp.Process(target = sim.start_run, args=(n,)))
#       
#    for P in Processes:
#        P.start()
#    for P in Processes:
#        P.join()
    
#    n = np.array([i+1 for i in range(len(velocities))])
#    sim.plot_pos(n)
#    sim.plot_p(n)
#    sim.plot_f(n)
###########################################