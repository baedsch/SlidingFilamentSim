import numpy as np
from numpy import random
import time as tme
from simulation3 import simulation#, Consumer#, varstep_sim

#initialization

#################################################################################################|

mode = 'springControl' #choose from ['vControl', 'fControl', ]
option = '' #for fControl choose from xy			
				#for vControl choose from 													
				#						-> poly: specify coefficients						
				#						-> step: specify n_elem, n_jumps, min_val, max_val		
name = 'spring_test_newInt_vRange__k_5__ns_5e4'													
																								
#store data in ram / write them to text files?													
s_store = False																					
p_store = True																					
f_store = False																					
sum_f_store = True																				
pos_store = True																			     	
writeText = False																			
																								
#most important parameters																		
n_heads = int(1e2)																			
n_steps = int(5e4)																				
d_t = 5e-3																						
bta = 2.																						
k = 10.
k_on = 10.																						
th = 0.01																						
t0 = 0.																						
d = 2.																				
random.seed(121155)																				
																								
#parameters for fControl																		    
loadF = [10. for i in range(3)]																
																								
#parameters for vControl																		    
	#-> poly option																				
v_Coeff = [[10., 0.] for i in range(3)] #example for constant velocity of 1.										
																								
	#-> step option																				
#step_n_jumps = int(n_steps / 1000)															    
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']												
#step_n_jumps = [25.,25]																		    
#step_min_val = [0., 0.05]																		
#step_max_val = [0.05, 0.2]
step_n_jumps = [25., 25.]																		    
step_min_val = [0., 0.05]																		
step_max_val = [0.05, 0.2]
																		
#step_factors = [1, 3]																		    
if not len(step_min_val) == len(step_max_val) and mode == 'vControl': 		                    
	raise ValueError("step_min_val, step_max_val and step_factors must have the same length")																							#|
           
#parameters for springControl
k_pull = 5
v_pull = 5

velocities = [0.5*i for i in range(5)]
                         																
#configure mulitprocessing																		
n_cores = 8																						
																								
																								
#################################################################################################|

########### MULTIPLE RUN, GIVEN SPRING VELOCITY SNIPPET

sim = simulation(mode,
    n_steps,
    n_heads,
    name = name,
    option = option,
    s_store = s_store,
    f_store = f_store,
    f_sum_store = sum_f_store,
    pos_store = pos_store,
    writeText = writeText,
    bta = bta,
    k = k,
    k_on = k_on,
    th = th,
    d_t = d_t,
    d = d,
    k_pull = k_pull,
    v_pull = v_pull)

for v in velocities:
	sim.add_run(v_pull=v)
n = np.array([i+1 for i in range(len(velocities))])
for rn in n:
	sim.start_run(rn)
#	sim.sum_up_P(rn)
sim.plot_pos(n)
sim.plot_p(n)
sim.plot_f(n)
###########################################