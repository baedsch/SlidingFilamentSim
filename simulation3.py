import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from numpy import random
import os
import platform
import time as tme
import time
#import threading as th
import multiprocessing as mp

#import tqdm
import bisect
import itertools
import collections
import helper3 as h




class simulation:
    #initialization:
    #kwargs: mode, seed, d, s, p, d_t, t0, k, bta, values_p, probabilities_p, vCoefficients (function of v polynomial)
    #checks, which variables should be stored, by default, everything is stored
    def __init__(self, mode, n_steps=int(1e4), n_heads=int(1e3), name="default", **kwargs):
        #these are the parameters
        #self.variables GET CHANGED WHITHIN FCN SCOPES!! no globaliyation needed (self. is like global)
        #most values are stored np.arrays nested in lists: for each run, append new array to list
        #within the arrays, axis 0 contains the time, axis 1 stands for the head index
        self.n_sim = [0]
        self.name = str(name)
        self.run = 0
        self.n_steps = [n_steps]
        self.n_heads = [n_heads]
        if mode in ['vControl', 'fControl', 'springControl']: self.mode = [mode]
        else: raise ValueError('Wrong mode chosen')
        self.option = [kwargs.get('option', '')]
        if not kwargs.get('option', '') in ['poly', 'step', 'const', '']: raise ValueError('Wrong option chosen')
        self.args_passed = [kwargs]

        self.d = [kwargs.get('d', 2.)]
        self.d_t = [kwargs.get('d_t', 5e-3)]
        self.seed = [kwargs.get('seed', False)]

        self.loadF = [kwargs.get('loadF', 20.)]
        self.bta = [kwargs.get('bta', 2.)]
        self.k = [kwargs.get('k', 10.)]
        self.k_on = [kwargs.get('k_on', 10.)]
        self.th = [kwargs.get('th', 0.0001)]
        self.n_div = [kwargs.get('n_div', 10)]
        self.k_pull = [kwargs.get('k_pull', 30)]
        self.v_pull = [kwargs.get('v_pull', 5)]
        self.min_max_k_min = [{}]
        self.min_max_k_plus = [{}]
        #values should look like this: [c_0, c_1, ..., c_n+1] for n^th polynomial
        self.v_Coeff = [kwargs.get('v_Coeff', [0])]
        self.v = [kwargs.get('v', 0)]
        self.repetitions = [kwargs.get('repetitions', 1)]
        self.step_n_jumps = [kwargs.get('step_n_jumps', int(n_steps / 1000))]
        self.step_min_val = [kwargs.get('step_min_val', 0.)]
        self.step_max_val = [kwargs.get('step_max_val', 100.)]

        self.values_p = [kwargs.get('values_p', [1, 0])]
        self.probabilities_p = [kwargs.get('probabilities_p', [0.7, 0.3])]
        self.s = [kwargs.get('s', h.random_cont(self.n_heads[self.run], -self.d[self.run]/2, self.d[self.run]/2))]
        self.p = [kwargs.get('p', h.random_discrete(self.n_heads[self.run], probabilities=self.probabilities_p[self.run]))]
        self.t0 = [kwargs.get('t0', 0.)]
        if mode in ['vControl']: self.t = [np.array([i*self.d_t[self.run] + self.t0[self.run] for i in range(self.n_steps[self.run])])]
        if mode in ['fControl', 'springControl']: self.t = [np.zeros(n_steps)]
#        self.sum_P = [np.array([])]
        self.n_neighbours = [h.find_neighbours(self.th[self.run], self.d[self.run], self.bta[self.run], self.k[self.run], self.k_on[self.run])]

        #declare storage of desired Variables out of s,p,f,sum_F
        self.store = {}
        self.store['s'] = kwargs.get('s_store',True)
        if self.store.get('s'): self.S = [np.zeros((n_steps, n_heads))]
        self.store['pos'] = kwargs.get('pos_store',True)
        if self.store.get('pos') or self.mode[self.run] in ['fControl', 'springControl']: self.Pos = [np.zeros((n_steps,1))]
        self.store['p'] = kwargs.get('p_store',True)
        if self.store.get('pos_pull') or self.mode[self.run] in ['springControl']: self.Pos_pull = [np.zeros((n_steps,1))]
        self.store['p'] = kwargs.get('p_store',True)
        if self.store.get('p'): self.P = [np.zeros((n_steps, n_heads))]
        self.store['f'] = kwargs.get('f_store',True)
        if self.store.get('f'): self.F = [np.zeros((n_steps, n_heads))]
        self.store['sum_f'] = kwargs.get('sum_f_store',True)
        if self.store.get('sum_f'): self.sum_F = [np.zeros((n_steps,1))]
        self.store['sum_p'] = kwargs.get('sum_p_store',True)
        if self.store.get('sum_p'): self.sum_P = [np.zeros((n_steps,1))]

        #write the stored variables into textfiles?
        self.writeText = kwargs.get('writeText', False)
        #make directory
        directory = tme.strftime("%Y_%m_%d__%H;%M;%S_", tme.gmtime()) + self.name
        try:
            if platform.system() == 'Windows':
                path = kwargs.get('path', 'C:/Users/baumgartner/sharePETER/sharePETER/res/' + directory)
                os.makedirs(path)
                os.chdir(path)
            if platform.system() == 'Linux':
                path = kwargs.get('path', '/res/' + directory)
                os.makedirs('.' + path)
                os.chdir(path)
        except:
           raise ValueError("Directory already existing, please choose unique simulation name! Otherwise data might get overwritten")
       
        with open('SimulationParameters', 'w') as file:
            file.write('{}:'.format(mode))
            file.write('Iteration Number per run: {}'.format(self.n_steps[0]))
            file.write('Number of trajectories at same velocity: {} (for getting statistical significance)'.format(self.repetitions[0]))
            file.write('kon {}'.format(self.k_on[0]))
            file.write('k0 {}'.format(self.k[0]))
            file.write('energyconstant (beta) {}'.format(self.bta[0]))
            file.write('d {}'.format(self.d[0]))
            file.write('neighbourcriterion {}'.format(self.th[0]))
            file.write('n heads {}'.format(self.n_heads[0]))
            file.write('initialpvector (pdetached, pattached) {}'.format(self.probabilities_p[0]))
            if mode == 'vControl': file.write('v prescribed from v={} to v={} (adapt the step such that it takes a reasonable time)'.format(self.step_min_val[0], self.step_max_val[0]))
            if mode == 'fControl': file.write('f prescribed from f={} to f={} such that you get total breakdown of all the bridges quickly in the simulation (adapt the step such that it takes a reasonable time)'.format(self.step_min_val[0], self.step_max_val[0]))
        
        #if filament breaks
        self.breakindex = [-1]
        
        #for step simulation in one run
        self.f_Sum_axis_added = {}
        self.v_axis_added = {}
        self.f_Sum_axis = [[]]
        self.v_axis = [[]]
        
        #PART FOR MULTIPLE SIMULATIONS -> add_simulation
        self.f_Sum_axis_runs = []
        self.v_axis_runs = []
    
    #similar to __init__
    def add_run(self, **kwargs):
        run = self.run
        self.n_sim.append(kwargs.get('n_sim', self.n_sim[self.run]))
        self.n_steps.append(kwargs.get('n_steps', self.n_steps[self.run]))
        self.n_heads.append(kwargs.get('n_heads', self.n_heads[self.run]))
        mode = kwargs.get('mode', self.mode[self.run])
        if mode in ['vControl', 'fControl', 'springControl']: self.mode.append(kwargs.get('mode', self.mode[self.run]))
        else: raise ValueError('Wrong mode chosen')
        self.option.append(kwargs.get('option', self.option[self.run]))
        if not self.option[-1] in ['poly', 'step', 'const', '']: raise ValueError('Wrong option chosen')
        self.args_passed.append(kwargs)

        self.d.append(kwargs.get('d', self.d[self.run]))
        self.d_t.append(kwargs.get('d_t', self.d_t[self.run]))
        self.seed.append(kwargs.get('seed', self.seed[self.run]))

        self.loadF.append(kwargs.get('loadF', self.loadF[self.run]))
        self.bta.append(kwargs.get('bta', self.bta[self.run]))
        self.k.append(kwargs.get('k', self.k[self.run]))
        self.k_on.append(kwargs.get('k_on', self.k_on[self.run]))
        self.th.append(kwargs.get('th', self.th[self.run]))
        self.n_div.append(kwargs.get('n_div', self.n_div[self.run]))
        self.k_pull.append(kwargs.get('k_pull', self.k_pull[self.run]))
        self.v_pull.append(kwargs.get('v_pull', self.v_pull[self.run]))
        self.min_max_k_min.append({})
        self.min_max_k_plus.append({})
        self.v_Coeff.append(kwargs.get('v_Coeff', self.v_Coeff[self.run]))
        self.v.append(kwargs.get('v', self.v[self.run]))
        self.repetitions.append(kwargs.get('repetitions', self.repetitions[self.run]))
        self.step_n_jumps.append(kwargs.get('step_n_jumps', int(self.n_steps[run] / 1000)))
        self.step_min_val.append(kwargs.get('step_min_val', 0.))
        self.step_max_val.append(kwargs.get('step_max_val', 100.))

        self.values_p.append(kwargs.get('values_p', self.values_p[self.run]))
        self.probabilities_p.append(kwargs.get('probabilities_p', self.probabilities_p[run]))
        self.t0.append(kwargs.get('t0', self.t0[self.run]))
        
        self.run += 1
        run = self.run
        self.s.append(kwargs.get('s', h.random_cont(self.n_heads[self.run], -self.d[self.run]/2, self.d[self.run]/2)))
        self.p.append(kwargs.get('p', h.random_discrete(self.n_heads[self.run], values=self.values_p[self.run], probabilities=self.probabilities_p[self.run])))
        

        if self.mode[self.run] in ['vControl']: self.t.append(np.array([i*self.d_t[self.run] + self.t0[self.run] for i in range(self.n_steps[self.run])]))
        elif self.mode[self.run] in ['fControl', 'springControl']: self.t.append(np.zeros(self.n_steps[run]))
#        self.sum_P.append(np.array([]))
        self.n_neighbours.append(h.find_neighbours(self.th[self.run], self.d[self.run], self.bta[self.run], self.k[self.run], self.k_on[run]))


        #append new slots for storage of desired Variables out of s,p,f,sum_F
        if self.store.get('s'): self.S.append(np.zeros((self.n_steps[run], self.n_heads[run])))
        if self.store.get('pos') or self.mode[self.run] in ['fControl', 'springControl']: self.Pos.append(np.zeros((self.n_steps[self.run],1)))
        if self.store.get('pos_pull') or self.mode[self.run] in ['springControl']: self.Pos_pull.append(np.zeros((self.n_steps[self.run],1)))
        if self.store.get('p'): self.P.append(np.zeros((self.n_steps[run], self.n_heads[run])))
        if self.store.get('f'): self.F.append(np.zeros((self.n_steps[run], self.n_heads[run])))
        if self.store.get('sum_f'): self.sum_F.append(np.zeros((self.n_steps[run],1)))
        if self.store.get('sum_p'): self.sum_P.append(np.zeros((self.n_steps[run],1)))

        self.breakindex.append(-1)
        self.f_Sum_axis.append([])
        self.v_axis.append([])
    
    #multiple simulations can be used to organize simulations, new simulation will overwrite all data but the listed down there
    #attention!!!
    #might overwright data
#    def add_simulation(self):
#        self.n_sim += 1
#        self.run = 0
#        #for step simulation in one-run-per-param
#        self.f_Sum_axis_runs.append([])
#        self.v_axis_runs.append([])
    
    # initialize a randomly chosen s (resp. p) as a startvector
    def init_s_rand(self):
        self.s[self.run] = h.random_cont(self.n_heads[self.run], -self.d[self.run]/2, self.d[self.run]/2)
        return self.s[self.run]
    def init_p_rand(self):
        self.p[self.run] = h.random_discrete(self.n_heads[self.run], self.values_p[self.run], self.probabilities_p[self.run])
        return self.p[self.run]
    
    #this is to get the upper and lower boundaries for acceptance rejection method
    def init_min_max_k(self, run):
        self.min_max_k_min[run] = h.min_max_min(self.n_div[run], self.bta[run], self.k[run])
        self.min_max_k_plus[run] = h.min_max_min(self.n_div[run], self.d[run], self.bta[run], self.k[run], self.k_on[run], self.n_neighbours[run])
        
    
    #################################################MAIN UPDATE METHODS
    #main update method for mode==vControl
    def update_vC(self, s, p, d, r01, run, step, bta, k):
        #update position according to option
        if self.option[run] == 'poly':
            v = np.polynomial.polynomial.polyval(step, self.v_Coeff[run])
        elif self.option[run] == 'step':
            v = h.step_fcn(self.n_steps[run], self.step_n_jumps[run], self.step_min_val[run], self.step_max_val[run], step)
        elif self.option[run] == 'const':
            v = self.v[run]
        else:
            raise ValueError("error! option: {}".format(self.option[run]))
        s = h.translate_def(s, v * self.d_t[run])

        #case: unbound
        if h.unitize(p) == 0:
            K_plus = h.k_plus(h.s_row(h.wrapping(s, self.d[run]), self.n_neighbours[run], self.d[run]), p, self.d[run], self.bta[run], self.k[run], self.k_on[run])
#            K_sum = h.k_plus_sum(s, p, self.d[run], self.bta[run], self.k[run], self.k_on[run], self.n_neighbours[run],  w=True)
            K_sum = sum(K_plus)
            if r01 <= K_sum * self.d_t[run]:
                p = h.det_p(h.p_row(self.n_neighbours[run]), K_plus, K_sum)
                s = h.wrapping(s, self.d[run])
                if p > 0: s += (p - 1) * d
                elif p < 0: s += p * d
                else: raise ValueError("something wrong with p")
        #case: bound
        else:
            if r01 <= h.k_min(h.force(s, p, self.d[run]), p, bta, k) * self.d_t[run]:
                p = 0

        return s, p
    #vectorize it in order to enhance performance
    updateV_vC = np.vectorize(update_vC)

    #main update method for mode==fControl
    #K_plus is the row of the k_plus_matrix corresponding to the updated head
    #K_sum is the sum of this row
    def update_fC(self, s, p, d, K_plus, K_sum, run):
        #case: unbound
        if h.unitize(p) == 0:
            p = h.det_p(h.p_row(self.n_neighbours[run]), K_plus, K_sum)
            #wrapping in here is mandatory for force calclation
            s = h.wrapping(s, self.d[run])
            if p > 0: s += (p - 1) * d
            elif p < 0: s += p * d
            else: raise ValueError("something wrong with p")
        
        #case: bound
        else: p = 0

        return s, p
    #vectorize it in order to enhance performance
    updateV_fC = np.vectorize(update_fC)
    
    def update_sC(self, s, p, d, K_plus, K_sum, run):
        #case: unbound
        if h.unitize(p) == 0:
            p = h.det_p(h.p_row(self.n_neighbours[run]), K_plus, K_sum)
            #wrapping in here is mandatory for force calclation
            s = h.wrapping(s, d)
            if p > 0: s += (p - 1) * d
            elif p < 0: s += p * d
            else: raise ValueError("something wrong with p")
        
        #case: bound
        else: p = 0

        return s, p
    
    #MAIN UPDATE METHODS#END############################################
    
    #run actual simulation
    def start_run(self, runs):
        #this loop enables running several runs with one function call -> Multiprocessing
        if isinstance(runs, int): runs = [runs]
        for run in runs:
                
            s = self.s[run]
            p = self.p[run]
            sum_F = 0.
            pos = 0.
            pos_pull = 0.
            t = 0.
    
            #~~~~~~~~~~~~~~~~SETUP INITIAL CONDITIONS BEFORE ITERATION~~~~~~~~~~~~~~~~~~~~~~
            #random numbers for udating (for 'fControl', we need two sets, for springControl three)
            if self.mode[run] in ['vControl']:
                print('rand creates')
                rand01 = random.random_sample((self.n_steps[run], self.n_heads[run]))
            
            if self.mode[run] in ['fControl']:
                rand011 = random.random_sample(self.n_steps[run])
                rand012 = random.random_sample(self.n_steps[run])
                n_dd = 0
                #stretching to desired force
                f = h.forceV(s, p, self.d[run])
    
                sum_F = sum(f)
                #~ print sum_F
                n_att = sum(h.unitizeV(p))
                #~ print n_att
                displ = -(sum_F - self.loadF[run]) / n_att
                #~ print displ
                s = h.translateV(s, displ)
                pos += displ
            
            if self.mode[run] in ['springControl']:
                n_k_p_u = 0
                n_k_p_nu = 0
                n_k_m_u = 0
                n_k_m_nu = 0
                n_dd = 0
                
                rand011 = random.random_sample(self.n_steps[run])
                rand012 = random.random_sample(self.n_steps[run])
                rand013 = random.random_sample(self.n_steps[run])
                k_min_max = h.get_max_k_min(self.bta[run], self.k[run])
                print('k_min_max:', k_min_max)
                k_plus_max = h.get_max_k_plus_sum(self.d[run], self.bta[run], self.k[run], self.k_on[run], self.n_neighbours[run])
                print('k_plus_max:', k_plus_max)
                
                #stretching to desired force
                f = h.forceV(s, p, self.d[run])
    
                sum_F = sum(f)
                print(sum_F)
                n_att = sum(h.unitizeV(p))
                print(n_att)
                displ = + sum_F / (n_att + self.k_pull[run])
                print(displ)
                s = h.translateV(s, displ)
                pos += displ
            #~~~~~~~~~~~~~~END SETUP INITIAL CONDITIONS BEFORE ITERATION~~~~~~~~~~~~~~~~~~~~~
            
            print('Run {}'.format(run))
            
            #loop over all iteration steps (install tqdm for progress bar)
        
#            for i in tqdm.tqdm(range(self.n_steps[run])):
            for i in range(self.n_steps[run]):

                #calculate force, f contains all forces of each head, sum_F is total load on filament
                f = h.forceV(s, p, self.d[run])
                sum_F = sum(f)
                sum_P = sum(h.unitize(p))
    
                #store Variables s, p
                if self.store.get('s'): self.S[run][i] = s
                if self.store.get('pos'): self.Pos[run][i] = pos
                if self.store.get('p'): self.P[run][i] = p
                if self.store.get('f'): self.F[run][i] = f
                if self.store.get('sum_f'): self.sum_F[run][i] = sum_F
                if self.store.get('sum_p'): self.sum_P[run][i] = sum_P
    
                #update s and p according to mode
                #self.s and self.p remain unchanged!
                ##begin case vControl
                if self.mode[run] == 'vControl':
                    s, p = self.updateV_vC(self, s, p, self.d[run], rand01[i], run, i, self.bta[run], self.k[run])
                    if self.option[run] == 'poly':
                        displ = np.polynomial.polynomial.polyval(i, self.v_Coeff[run])
                    elif self.option[run] == 'step':
                        displ = h.step_fcn(self.n_steps[run], self.step_n_jumps[run], self.step_min_val[run], self.step_max_val[run], i)
                    elif self.option[run] == 'const':
                        displ = self.v[run]
                    pos += displ *  self.d_t[run]
                ##end case vControl
                
                ##begin case fControl
                elif self.mode[run] == 'fControl':
                    #immediate detachment if heads too far away from binding site
                    for hi in range(self.n_heads[run]):
                        if h.unitize(p[hi]) and abs(s[hi]) > (1 + self.k[run]) / 2:
                             p[hi] = 0
                             n_dd += 1
                    
                    #calculate number of attached heads BEFORE AND AFTER UPDATE STEP
                    n_att = sum(h.unitizeV(p))
                    if n_att == 0:
                        print("connection broke! This happened at iterationstep {}".format(i))
                        self.t[run][i:] = t
                        self.breakindex[run] = i
                        break
                    
                    #here, probabilities for attaching, detaching and the corresponding wating time tau are calculated
                    k_plus_sum = h.k_plus_sum(s, p, self.d[run], self.bta[run], self.k[run], self.k_on[run], self.n_neighbours[run],  w=True)
                    k_min = h.k_minV(s, p, self.bta[run], self.k[run])
                    k = k_plus_sum + k_min
                    k_a = list(itertools.accumulate(k))
                    k_sum = sum(k_plus_sum + k_min)
    
                    #find the result for the waiting time
                    tau = -1 / k_sum * np.log(rand011[i])
                    
                    #get index of selected head
                    min_index = bisect.bisect_left(k_a, k_sum * rand012[i])
                    s_row = h.s_row(h.wrapping(s[min_index],self.d[run]), self.n_neighbours[run], self.d[run])
                    k_plus_row = h.k_plus_matrix(s_row, p[min_index], self.d[run], self.bta[run], self.k[run], self.k_on[run])
                    
                    #update this head, calculate the force difference
                    s_i, p_i = s[min_index], p[min_index]
                    
                    #in case of binding, s gets wrapped and positioned acc to p
                    s_upd, p_upd = self.update_fC(s_i, p_i, self.d[run], k_plus_row, sum(k_plus_row), run)
                    s[min_index] = s_upd
                    p[min_index] = p_upd
                    
                    #calculate force vector of updated state and the delta
                    f = h.forceV(s, p, self.d[run])
                    f_delta = sum(f) - self.loadF[run]
                    
                    #calculate corresponding displacement
                    displ = -f_delta / float(n_att)
    
                    #updating filament position and s positions
                    pos += displ
                    s = h.translateV(s, displ)
    
                    #updating the elapsed time
                    #~ print tau_min
                    t += tau
                    self.t[run][i] = t
                    
                    #calculate number of attached heads BEFORE AND AFTER UPDATE STEP
                    n_att = sum(h.unitizeV(p))
                    
                    #check again, if connection broke
                    if n_att == 0:
                        print("connection broke! This happened at iterationstep {}".format(i))
                        self.t[run][i:] = t
                        self.breakindex[run] = i
                        break
                ##end case fControl
                
                ##begin case springControl
                elif self.mode[run] == 'springControl':
                    #immediate detachment if heads too far away from binding site
                    for hi in range(self.n_heads[run]):
                        if h.unitize(p[hi]) and abs(s[hi]) > (1 + self.k[run]) / 2:
                             p[hi] = 0
                             n_dd += 1
                    
                    #calculate number of attached heads BEFORE AND AFTER UPDATE STEP
                    n_att = sum(h.unitizeV(p))
                    if n_att == 0:
                        print("connection broke! This happened at iterationstep {}".format(i))
                        self.t[run][i:] = t
                        self.breakindex[run] = i
                        break
                    
                    #upper prospensity bound
                    k_upper = np.zeros(self.n_heads[run])
                    for hi in range(self.n_heads[run]):
                        #case: unbound
                        if h.unitize(p[hi]) == 0:
                            k_upper[hi] = k_plus_max
                        #case: bound
                        if h.unitize(p[hi]) == 1:
                            k_upper[hi] = k_min_max
                    k_upper_sum = sum(k_upper)
                    
                    #calc tau
                    tau = - np.log(rand011[i]) / k_upper_sum
                    t += tau
                    self.t[run][i] = t
                    
                    ##displacements during waiting time
                    pos_pull_upd = pos_pull + self.v_pull[run] * tau
                    n_att = sum(h.unitizeV(p))
                    
                    #attention! could be too many subtractions
                    delta_pos = (self.k_pull[run] * pos_pull_upd - sum_F + n_att * pos) / (self.k_pull[run] + n_att) - pos
                    #apply translation
                    s = h.translateV(s, delta_pos)
                    pos += delta_pos
                    pos_pull = pos_pull_upd
                    
                    #get index of update candidate
                    k_upper_accumulated = np.add.accumulate(k_upper)
                    index = bisect.bisect_left(k_upper_accumulated, k_upper_sum * rand012[i])
                    n_att = sum(h.unitizeV(p))
                    
                    s_i, p_i = s[index], p[index]
                    
                    #verify the head selection: if RN < prob, head is accepted, prob is integrated rate k(t) dt
                    ##begin case: distinguish between
                        # bound -> unbound case (simple)
                        # unbound -> bound case (complicated, because sum over complicated integrals)
                    #incrementing variables e.g. n_k_p_u for debugging
                    #case: unbound
                    if h.unitize(p[index]) == 0:
                        ##calculate probability ratio in order to verify head selection
                        v = delta_pos / tau
                        #attention! could be too many subtractions (optimization could be good...)
                        prob = h.int_k_plus_sum(s[index] - delta_pos, s[index], self.d[run], self.bta[run], self.k[run], self.k_on[run], self.n_neighbours[run], v)
                        if rand013[i] < prob / (k_plus_max * tau):
                            n_k_p_u +=1
                            #rows are line vector containing head pos relative to neighbours taken into account and derived values
                            s_row = h.s_row(h.wrapping(s[index],self.d[run]), self.n_neighbours[run], self.d[run])
                            k_plus_row = h.k_plus_matrix(s_row, p[index], self.d[run], self.bta[run], self.k[run], self.k_on[run])
                            #update step
                            s[index], p[index] = self.update_sC(s[index], p[index], self.d[run], k_plus_row, sum(k_plus_row), run)
                        else: n_k_p_nu +=1
                    
                    #case: bound
                    else:
                        ##calculate probability ratio in order to verify head selection
                        #average velocity during waiting time(for integration in prob)
                        v = delta_pos / tau
                        prob = h.int_k_min(s[index] - delta_pos, s[index], p[index], self.bta[run], self.k[run], v)
                        prob = abs(prob)
                    
                        if rand013[i] < abs(prob) / (k_min_max * tau):
                            n_k_m_u +=1
                            #update step
                            s[index], p[index] = self.update_sC(s[index], p[index], self.d[run], 0, 0, run)
                        else: n_k_m_nu +=1
                    ##end case
                    
                    s_upd, p_upd = s[index], p[index]
                    n_att_upd = sum(h.unitizeV(p))
                    
                    #jump of filament due to attachment / detatchment
                    f = h.forceV(s, p, self.d[run])
                    sum_F = sum(f)
                    pos_upd = (self.k_pull[run] * pos_pull - sum_F + n_att_upd * pos) / (n_att_upd + self.k_pull[run])
                    jump = pos_upd - pos
                    #apply translation (jump)
                    s = h.translateV(s, jump)
                    pos = pos_upd
                    
                    #check again, if connection broke
                    n_att = sum(h.unitizeV(p))
                    if n_att == 0:
                        print("connection broke! This happened at iterationstep {}".format(i))
                        self.t[run][i:] = t
                        self.breakindex[run] = i
                        break
                ##end case springControl
                
            #option for debugging / understand what was going on during simulation
            if self.mode[run] == 'springControl' and False:        
                print('n_k_p_u:', n_k_p_u)
                print('n_k_p_nu:', n_k_p_nu)
                print('n_k_m_u:', n_k_m_u)
                print('n_k_m_nu:', n_k_m_nu)          
            if self.mode[run] in ['springControl', 'fControl']:print('n_dd:', n_dd)          
            
            #preparing the time vector to be added in first column of each file
            t_w = self.t[run][np.newaxis]
            t_w = t_w.T
            
            if self.option[run] in ['const', 'poly'] and self.mode[run] == 'vControl':
                self.average_norm_force_single(run)
            if self.mode[run] == 'fControl':
                self.average_velocity_single(run)
            
            #writing the results to textfiles
            if self.writeText:
                if self.mode[run] == 'vControl' and self.option[run] == 'const':
                        np.savetxt('time_{}prescribed_{}_Run{}.dat'.format(self.mode[run][0], self.v[run], self.n_sim[run]), t_w, header='time')
                elif self.mode[run] == 'fControl':
                    np.savetxt('time_{}prescribed_{}_Run{}.dat'.format(self.mode[run][0], self.loadF[run], self.n_sim[run]), t_w, header='time')
                else:
                    np.savetxt('time_{}_{}.dat'.format(run, self.n_sim[run]), t_w, header='time')
                    
                if self.store.get('s'):
                    S_out = np.concatenate((t_w, self.S[run]), axis=1)
                    np.savetxt('S_{}_{}.dat'.format(run, self.n_sim[run]), S_out, header='time, Position of each head')
                if self.store.get('pos'):
                    Pos_out = np.concatenate((t_w, self.Pos[run]), axis=1)
                    np.savetxt('Pos_{}_{}.dat'.format(run, self.n_sim[run]), Pos_out, header='time, Position of filament')
                if self.store.get('p'):
                    if self.mode[run] == 'vControl' and self.option[run] == 'const':
                        np.savetxt('p_{}prescribed_{}_Run{}.dat'.format(self.mode[run][0], self.v[run], self.n_sim[run]), self.P[run], header='time, force applied by each head')
                    elif self.mode[run] == 'fControl':
                        np.savetxt('p_{}prescribed_{}_Run{}.dat'.format(self.mode[run][0], self.loadF[run], self.n_sim[run]), self.P[run], header='time, force applied by each head')
                    else:
                        P_out = np.concatenate((t_w, self.P[run]), axis=1)
                        np.savetxt('P_{}_{}.dat'.format(run, self.n_sim[run]), P_out, header='time, binding state of each head')
                if self.store.get('f'):
                    if self.mode[run] == 'vControl' and self.option[run] == 'const':
                        np.savetxt('f_{}prescribed_{}_Run{}.dat'.format(self.mode[run][0], self.v[run], self.n_sim[run]), self.F[run], header='time, force applied by each head')
                    elif self.mode[run] == 'fControl':
                        np.savetxt('f_{}prescribed_{}_Run{}.dat'.format(self.mode[run][0], self.loadF[run], self.n_sim[run]), self.F[run], header='time, force applied by each head')
                    else:
                        F_out = np.concatenate((t_w, self.F[run]), axis=1)
                        np.savetxt('sum_F_{}_Run_{}.dat'.format(run, self.n_sim[run]), F_out, header='time, force applied by each head')
                if self.store.get('sum_f'):
                    if self.mode[run] == 'vControl' and self.option[run] == 'const':
                        np.savetxt('ftotal_{}prescribed_{}_Run{}.dat'.format(self.mode[run][0], self.v[run], self.n_sim[run]), self.sum_F[run], header='time, force applied by each head')
                    elif self.mode[run] == 'fControl':
                        np.savetxt('ftotal_{}prescribed_{}_Run{}.dat'.format(self.mode[run][0], self.loadF[run], self.n_sim[run]), self.sum_F[run], header='time, force applied by each head')
                    else:
                        sum_F_out = np.concatenate((t_w, self.sum_F[run]), axis=1)
                        np.savetxt('sum_F_{}_{}.dat'.format(run, self.n_sim[run]), sum_F_out, header='time, force applied by each head')
                if self.store.get('sum_p'):
                    if self.mode[run] == 'vControl' and self.option[run] == 'const':
                        np.savetxt('ptotal_{}prescribed_{}_Run{}.dat'.format(self.mode[run][0], self.v[run], self.n_sim[run]), self.sum_P[run], header='heads attached')
                    elif self.mode[run] == 'fControl':
                        np.savetxt('ptotal_{}prescribed_{}_Run{}.dat'.format(self.mode[run][0], self.loadF[run], self.n_sim[run]), self.sum_P[run], header='heads attached')
                    else:
                        sum_F_out = np.concatenate((t_w, self.sum_F[run]), axis=1)
                        np.savetxt('sum_P_{}_{}.dat'.format(run, self.n_sim[run]), sum_F_out, header='heads attached')
                
                #write simulation parameters
                if self.mode[run] == 'vControl' and self.option[run] == 'const':
                    with open('ParameterValues_{}_Run{}'.format(self.v[run], self.v[run])) as f:
                        f.write('Nrealisation kon timeconstant k0 energyconstant d0 Δt0 neighbourcriterion n initialpvector (pdetached, pattached) vprescribed\n')
                        f.write('{} {} {} {} {}	{}	{} {} {}	{} {}	{}\n'.format(self.n_steps[run], self.k_on[run], 'timeconstant', self.k[run], self.bta[run], self.d[run], self.d_t[run], self.th[run], 'n', 'initPVec', self.probabilities_p[run], self.v[run]))
                elif self.mode[run] == 'fControl':
                    with open('ParameterValues_{}_Run{}'.format(self.loadF[run], self.n_sim[run]), 'w') as f:
                        f.write('Nrealisation kon timeconstant k0 energyconstant d0 Δt0 neighbourcriterion n initialpvector (pdetached, pattached) fprescribed\n')
                        f.write('{} {} {} {} {}	{}	{} {} {}	{} {}	{}\n'.format(self.n_steps[run], self.k_on[run], 'timeconstant', self.k[run], self.bta[run], self.d[run], self.d_t[run], self.th[run], 'n', 'initPVec', self.probabilities_p[run], self.loadF[run]))
                else:
                    with open('ParameterValues_Run{}'.format(run)) as f:
                   
                        f.write('Nrealisation kon timeconstant k0 energyconstant d0 Δt0 neighbourcriterion n initialpvector (pdetached, pattached)\n')
                        f.write('{} {} {} {} {}	{}	{} {} {}	{} {}\n'.format(self.n_steps[run], self.k_on[run], 'timeconstant', self.k[run], self.bta[run], self.d[run], self.d_t[run], self.th[run], 'n', 'initPVec', self.probabilities_p[run]))
                

        return 1

    #takes the P array and returns the sum of bound heads for each timestep (after simulation)
    def sum_up_P(self, run):
        for p in self.P[run]:
            s = sum(h.unitize(p))
            self.sum_P[run] = np.append(self.sum_P[run], s)
        return self.sum_P[run]


    def plot_p(self, run, leg=False):
        if isinstance(run, int):
            run = [run]
            print('input into list')

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_title('Time evolution of bound heads ({} run(s)), mode: {}'.format(len(run), self.mode[run[0]]))
        ax.set_xlabel('time [s]')
        ax.set_ylabel('bound heads')

        for r in run:
            if len(self.args_passed[r]) == 1:
                lab = list(self.args_passed[r])[0]
                val = self.args_passed[r][lab]
#                print('t', max(self.t[r]))
#                print('p', self.sum_P[r])
                ax.plot(self.t[r], self.sum_P[r], linewidth=1.0, linestyle="-", label='{}={}'.format(lab, val))
            else: ax.plot(self.t[r], self.sum_P[r], linewidth=1.0, linestyle="-")
        if leg: ax.legend()
        plt.savefig('Sum_p.png', dpi=200)

    def plot_f(self, run, leg=False):
        if isinstance(run, int):
            run = [run]
            print('input into list')

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_title('Total load on filament ({} run(s)), mode: {}'.format(len(run), self.mode[run[0]]))
        ax.set_xlabel('time [s]')
        ax.set_ylabel('Total load')

        for r in run:
            if len(self.args_passed[r]) == 1:
                lab = list(self.args_passed[r])[0]
                val = self.args_passed[r][lab]
                ax.plot(self.t[r], self.sum_F[r], linewidth=1.0, linestyle="-", label='{}={}'.format(lab, val))
            else: ax.plot(self.t[r], self.sum_F[r], linewidth=1.0, linestyle="-")
        if leg: ax.legend()
        plt.savefig('Sum_f.png', dpi=200)

    def plot_f_norm(self, run, leg=False):
        if isinstance(run, int):
            run = [run]
            print('input into list')

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_title('Total load on filament ({} run(s)), mode: {}'.format(len(run), self.mode[run[0]]))
        ax.set_xlabel('time [s]')
        ax.set_ylabel('Total load normalized')

        for r in run:
            if len(self.args_passed[r]) == 1:
                lab = list(self.args_passed[r])[0]
                val = self.args_passed[r][lab]
                ax.plot(self.t[r], self.sum_F[r] / float(self.n_heads[r]), linewidth=1.0, linestyle="-", label='{}={}'.format(lab, val))
            else: ax.plot(self.t[r], self.sum_F[r] / float(self.n_heads[r]), linewidth=1.0, linestyle="-")
        if leg: ax.legend()
        plt.savefig('Sum_f.png', dpi=200)
    
    #
    def plot_f_norm__v(self,leg=False, c='b'):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_title('Total load on filament (mode: {}'.format(self.mode[0]))
        ax.set_xlabel('velocity [m/s]')
        ax.set_ylabel('Total load normalized')
        

        color=c
        ax.scatter(self.v_axis_runs, self.f_Sum_axis_runs, color=color, s=5)
        
        plt.savefig('f_norm__v.png', dpi=200)
        print(self.f_Sum_axis_runs)
    
    #TBA normalization of force!!
    def plot_v__f_norm(self, leg=False, c='b'):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_title('Total load on filament (mode: {}'.format(self.mode[0]))
        ax.set_xlabel('Total load')
        ax.set_ylabel('velocity [m/s]')
        #ax.set_yscale('log')
        color = c
        
        ax.scatter(self.f_Sum_axis_runs, self.v_axis_runs, color=color, s=5)
        
        plt.savefig('v__f_norm.png', dpi=200)
    
    def plot_pos(self, run, leg=False):
        if isinstance(run, int):
            run = [run]
            print('input into list')

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_title('Position filament ({} run(s)), mode: {}'.format(len(run), self.mode[run[0]]))
        ax.set_xlabel('time [s]')
        ax.set_ylabel('position')

        for r in run:
            if len(self.args_passed[r]) == 1:
                lab = list(self.args_passed[r])[0]
                val = self.args_passed[r][lab]
                ax.plot(self.t[r], self.Pos[r], linewidth=1.0, linestyle="-", label='{}={}'.format(lab, val))
            else: ax.plot(self.t[r], self.Pos[r], linewidth=1.0, linestyle="-")
        if leg: ax.legend()
        plt.savefig('Pos.png', dpi=200)

    def plot_f_v_step(self, run, leg=False, legvar='', legval=[], c='b'):
        if isinstance(run, int):
            run = [run]
            print('input into list')

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_title('Load per head on filament ({} run(s)), mode: {}'.format(len(run), self.mode[run[0]]))
        ax.set_xlabel('velocity')
        ax.set_ylabel('Total load')
        val = 0

        for r in run:
            if not isinstance(c, str): color = c[r]
            else: color = c
            lab = legvar
            if leg: val = legval[run.index(r)]
            ax.scatter(self.v_axis_added[r], self.f_Sum_axis_added[r], color=color, s=5, label='{}={}'.format(lab, val))
        if leg: ax.legend()
        plt.savefig('Sum_f_v.png', dpi=200)

    def plot_v_step(self, run, leg=False, c='b'):
        if isinstance(run, int):
            run = [run]
            print('input into list')

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_title('Velocity of filament ({} run(s)), mode: {}'.format(len(run), self.mode[run[0]]))
        ax.set_xlabel('jump [arb. time]')
        ax.set_ylabel('velocity')

        for r in run:
            if len(self.args_passed[r]) == 1:
                lab = list(self.args_passed[r])[0]
                val = self.args_passed[r][lab]
                ax.plot(self.v_axis[r], color=c, label='{}={}'.format(lab, val))
            else: ax.plot(self.v_axis[r], color=c)
        if leg: ax.legend()
        plt.savefig('v.png', dpi=200)

    

    def create_threads(self, runs):
        self.jobs = [th.Thread(target=self.start_run, args=(i,)) for i in runs]
    def start_threads(self):
        for j in self.jobs:
            j.start()
        for j in self.jobs:
            j.join()
    
    #get the asymptotic force value for one run
    def average_norm_force_single(self, run, equilib_wait_frac=0.25):
        n = len(self.sum_F[run])
        print(sum(self.sum_F[run]))
        f_mean = np.mean(self.sum_F[run][(int(n * equilib_wait_frac)):])
        print(f_mean)
        self.f_Sum_axis_runs.append(f_mean / self.n_heads[run])
        self.v_axis_runs.append(self.v[run])
        
        while True: 
            try:
                with open('F_v.dat', 'a') as f:
                    f.write("{} {}\n".format(self.v[run], f_mean / self.n_heads[run]))
                break
            except: print("have to wait for textfile")
            
        return f_mean / self.n_heads[run], self.v[run]
    
   
    #get the asymptotic velocity value for one run
    def average_velocity_single(self, run, equilib_wait_frac=0.25):
        #case no break
        if self.breakindex[run] == -1: n = self.n_steps[run]
        #case break: only take data up to breakindex -2 into account
        else: n = self.breakindex[run] - 2#check?????????????????????????????????????????????????
        
        t = np.array(self.t[run][int(n * equilib_wait_frac):n])
        pos = self.Pos[run][int(n * equilib_wait_frac):n,0]
        
        popt, pcov = curve_fit(h.lin_fit, t, pos)
        v = popt[0]
        
        self.f_Sum_axis_runs.append(self.loadF[run])
        self.v_axis_runs.append(v)
        
        while True: 
            try:
                with open('v_F.dat', 'a') as f:
                    f.write("{} {}\n".format(self.loadF[run], v))
                break
            except: print("have to wait for textfile")
            
        return 1
    
	#get the asymptotic value for each step of v_step
    def average_norm_force(self, run, f_Sum, v_step_list, equilib_wait_frac=0.25):

        counter = collections.Counter(v_step_list)
        print(len(v_step_list))
        print(counter)
        v_axis = sorted([i for i in counter])
        f_Sum_axis = []
        index = 0
        for e in v_axis:
            n = counter.get(e)
            f = np.mean(f_Sum[(index + int(n * equilib_wait_frac)):(index + n)])
            f_Sum_axis.append(f / self.n_heads[run])
            index += n
        self.v_axis[run] = v_axis
        self.f_Sum_axis[run] = f_Sum_axis
        return v_axis, f_Sum_axis

    def join_steps(self, runs, index):
        runs.sort()

        self.f_Sum_axis_added[index] = []
        self.v_axis_added[index] = []

        for e in runs:
            self.f_Sum_axis_added[index] += self.f_Sum_axis[e]
            self.v_axis_added[index] += self.v_axis[e]
