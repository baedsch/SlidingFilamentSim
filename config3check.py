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
                    s=[-0.5093321181141888,
 -0.681112344058524,
 -0.43202495465981317,
 0.18293457146938286,
 -0.07648946416162383,
 0.0020131366483209234,
 0.7841875980533466,
 -0.836000543754047,
 0.8862758367210133,
 0.955243017716624,
 0.30520312457143906,
 -0.7085465721847561,
 0.10941748646789184,
 -0.8946966128222149,
 0.3431225721368336,
 -0.6103248804183066,
 -0.6960484816709815,
 0.05167743248455414,
 -0.805766551928107,
 -0.994261895312988,
 0.8344591021204155,
 0.26579268456393734,
 0.7500856423295263,
 0.17902426065431376,
 0.049952543761914825,
 -0.7399771824461432,
 0.8911813859325814,
 -0.3670218949323254,
 -0.5223444061111566,
 0.34433510589547867,
 0.3414039706038121,
 0.11184830003806967,
 0.4053260408003907,
 -0.21557540261892472,
 -0.1915320482997407,
 0.6610736513760944,
 -0.8373102647544042,
 0.0026973100555971463,
 -0.7816199229825243,
 0.9225985630042488,
 0.1527286424585581,
 -0.5054191388726386,
 0.8894578129744648,
 0.8592876228256374,
 -0.12785791063496355,
 0.9589658606470359,
 0.013603278692741139,
 -0.4806013853679625,
 0.07005950385046589,
 -0.2444753727770168,
 0.6084573043038168,
 -0.05073624008745581,
 -0.22270332344843347,
 -0.35624827117441304,
 0.15715937867216034,
 -0.24048648185094823,
 -0.040733292041214675,
 0.44285470846251185,
 -0.773277180886895,
 0.1047267059846444,
 -0.06091804486023378,
 0.5548010367895637,
 -0.7032465006254927,
 0.9327979001490219,
 0.034043987636752426,
 -0.40507341944797703,
 -0.8247322717155448,
 -0.3087786845253033,
 0.9102394086303232,
 0.7167205890608286,
 0.17575002158560338,
 0.3535255894105571,
 -0.4152856822858988,
 -0.69980763723219,
 -0.27043155572008004,
 0.017456030697633462,
 0.6194463131852475,
 -0.4137566355131126,
 -0.4441269188495873,
 0.0022962809325641764,
 0.04099026325252364,
 -0.5336165459618716,
 0.6085306884075035,
 0.5674714895429434,
 -0.8600798550584341,
 -0.23681256478990242,
 -0.3201486237050797,
 0.8976762901719453,
 -0.33937340483300105,
 0.19689259978727147,
 -0.495425478757904,
 -0.577326082305456,
 -0.6660334311243008,
 -0.4861756887615396,
 0.6758157766871551,
 0.1356911249572419,
 -0.7569920663189604,
 0.9943903827899554,
 -0.4206005447393233,
 0.5878676243794425],
                    p=[1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 0.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 0.0,
 1.0,
 0.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 0.0,
 0.0,
 1.0,
 0.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 0.0,
 1.0,
 0.0,
 1.0,
 0.0,
 1.0,
 0.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 0.0,
 1.0,
 0.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 0.0,
 1.0,
 0.0,
 1.0,
 0.0,
 1.0,
 0.0,
 1.0,
 1.0,
 0.0,
 1.0,
 1.0,
 1.0,
 0.0,
 1.0,
 0.0,
 1.0,
 1.0,
 1.0,
 0.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 1.0,
 0.0,
 1.0,
 1.0,
 0.0,
 0.0,
 1.0,
 0.0,
 1.0,
 1.0,
 1.0,
 0.0,
 0.0],
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