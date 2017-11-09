import numpy as np
from numpy import random
import time as tme
from simulation3 import simulation#, Consumer#, varstep_sim

#initialization

#################################################################################################|

mode = 'fControl' #choose from ['vControl', 'fControl', ]
option = 'poly' #for fControl choose from xy			
				#for vControl choose from 													
				#						-> poly: specify coefficients						
				#						-> step: specify n_elem, n_jumps, min_val, max_val		
name = 'matematica_crosscheck_f=10'													
																								
#store data in ram / write them to text files?													
s_store = False																					
p_store = True																					
f_store = False																					
sum_f_store = True																				
pos_store = True																			     	
writeText = True																				
																								
#most important parameters																		
n_heads = int(1e2)																			
n_steps = int(1e4)																				
d_t = 5e-3																						
bta = 2.																						
k = 10.
k_on = 10.																						
th = 0.01																						
t0 = 0.																						
d = 2.																				
random.seed(121155)																				
																								
#parameters for fControl																		    
loadF = [10. for i in range(2)]																
																								
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
                                    																
#configure mulitprocessing																		
n_cores = 8																						
																								
																								
#################################################################################################|

time_before = tme.time()

#sim = simulation(	mode,
#					n_steps,
#					n_heads,
#					name = name,
#					option = option,
#					loadF = loadF,
#					v_Coeff = v_Coeff,
#					step_n_jumps = step_n_jumps[0],
#					step_min_val = step_min_val[0],
#					step_max_val = step_max_val[0],
#					s_store = s_store,
#					f_store = f_store,
#					f_sum_store = sum_f_store,
#					pos_store = pos_store,
#					writeText = writeText,
#					bta = bta,
#					k = k,
#					k_on = k_on,
#					th = th,
#					d_t = d_t,
#					d = d[0])


########### trying mult. possibilities for initial p vectors
#probabilities_p = [[a / 10., round(1.- a / 10., 2)] for a in range(11)]
#for p in probabilities_p:
	#sim.add_run(probabilities_p=p)
############################################################



########### MULTIPLE RUN, DIFF LOAD SNIPPET
#for f in loadF:
	#sim.add_run(loadF=f)
#n = np.array([i+1 for i in range(len(loadF))])
#for rn in n:
	#sim.start_run(rn)
	#sim.sum_up_P(rn)
#sim.plot_pos(n)
#sim.plot_p(n)
#sim.plot_f(n)
###########################################

########### MULTIPLE RUN, GIVEN LOAD SNIPPET

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
					d = d)
for f in loadF:
	sim.add_run(loadF=f)
n = np.array([i+1 for i in range(len(loadF))])
for rn in n:
	sim.start_run(rn)
	sim.sum_up_P(rn)
sim.plot_pos(n)
sim.plot_p(n)
sim.plot_f(n)
###########################################

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
#                        k_on = k_on,
#					th = th,
#					d_t = d_t,
#					d = d)
#for f in v_Coeff:
#	sim.add_run(v_Coeff=f)
#n = np.array([i+1 for i in range(len(v_Coeff))])
#for rn in n:
#	sim.start_run(rn)
#	sim.sum_up_P(rn)
#sim.plot_pos(n)
#sim.plot_p(n)
#sim.plot_f(n)
###########################################

########### ONE SINGLE RUN SNIPPET
#sim.start_run(0)
#sim.sum_up_P(0)
#sim.plot_f(0)
#sim.plot_p(0)
#sim.plot_pos(0)

##################################

########### ONE SINGLE RUN SNIPPET WITH F-v fcn
#sim.start_run(0)
#sim.sum_up_P(0)
#sim.average_norm_force(0, sim.sum_F[0], h.step_list(n_steps, step_n_jumps[0], step_min_val[0], step_max_val[0])) #returns v, f axis, needs processing
#sim.plot_f(0)
#sim.plot_p(0)
#sim.plot_pos(0)
#sim.plot_f_v_step(0)
#sim.plot_v_step(0)
###############################################

########### MULTIPLE RUN SNIPPET WITH F-v fcn
#sim = simulation(	mode,
#					n_steps,
#					n_heads[0],
#					name = name,
#					option = option,
#					s_store = s_store,
#					f_store = f_store,
#					f_sum_store = sum_f_store,
#					pos_store = pos_store,
#					writeText = writeText,
#					bta = bta,
#					k = k,
#					th = th,
#					d_t = d_t,
#					d = d)
#z=1
#y=1
#indices = []
#for var1 in n_heads:
#    runs = []
#
#    for varis in range(len(step_max_val)):
#        sim.add_run(    n_heads = var1,
#                step_n_jumps = step_n_jumps[varis],
#                        step_min_val = step_min_val[varis],
#                        step_max_val = step_max_val[varis],
#                        )
#        sim.init_p_rand()
#        sim.init_s_rand()
#        sim.start_run(y)
#        sim.sum_up_P(y)
#
#        sim.average_norm_force(y, sim.sum_F[y], h.step_list(n_steps, step_n_jumps[varis], step_min_val[varis], step_max_val[varis]), equilib_wait_frac=0.25) #returns v, f axis, needs processing
#        runs.append(y)
#        y += 1
#
#    sim.join_steps(runs, z)
#
#    i=[v +z for v in range(len(step_max_val))]
#
#
#
#    indices.append(z)
#    z += 1
#sim.plot_f_v_step(indices, c=colors, leg=True, legvar='n_heads', legval=n_heads)
#print(sim.f_Sum_axis_added, sim.v_axis_added)
###############################################


time_after = tme.time()
print(time_after - time_before)
print('Seconds')


#if __name__ == '__main__':
#    # Establish communication queues
#    tasks = multiprocessing.JoinableQueue()
#    results = multiprocessing.Queue()
#
#    # Start consumers
#    num_consumers = multiprocessing.cpu_count() - 1
#    print('Creating {} consumers'.format(num_consumers))
#    consumers = [Consumer(tasks, results) for i in range(num_consumers)]
#    for w in consumers:
#        w.start()
#
#    # Enqueue jobs
#    num_jobs = len(v_Coeff)
#    n = np.array([i+1 for i in range(num_jobs)])
#    for r in n:
#        tasks.put(sim.run(r))
#
#    # Add a poison pill for each consumer
#    for i in range(num_consumers):
#        tasks.put(None)
#
#    # Wait for all of the tasks to finish
#    tasks.join()
#
#    # Start printing results
#    while num_jobs:
#        result = results.get()
#        print('Result:', result)
#        num_jobs -= 1