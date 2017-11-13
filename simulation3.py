import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
from numpy import random
#~ import random
import math
import os
import platform
import time as tme
import time
#import threading as th
import multiprocessing as mp
import multiprocessing
import tqdm
import bisect
import itertools
import collections
# in order to be able to use map with KeyWordArgs
from functools import partial
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
        self.step_n_jumps = [kwargs.get('step_n_jumps', int(n_steps / 1000))]
        self.step_min_val = [kwargs.get('step_min_val', 0.)]
        self.step_max_val = [kwargs.get('step_max_val', 100.)]

        self.values_p = [kwargs.get('values_p', [0., 1.])]
        self.probabilities_p = [kwargs.get('probabilities_p', [0.3, 0.7])]
        self.s = [kwargs.get('s', h.random_cont(self.n_heads[self.run], -self.d[self.run]/2, self.d[self.run]/2))]
        self.p = [kwargs.get('p', h.random_discrete(self.n_heads[self.run], probabilities=self.probabilities_p[self.run]))]
        self.t0 = [kwargs.get('t0', 0.)]
        if mode in ['vControl']: self.t = [np.array([i*self.d_t[self.run] + self.t0[self.run] for i in range(self.n_steps[self.run])])]
        if mode in ['fControl']: self.t = [np.zeros(n_steps)]
        self.sum_P = [np.array([])]

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

        self.f_Sum_axis_added = {}
        self.v_axis_added = {}
        self.f_Sum_axis = [[]]
        self.v_axis = [[]]
   
    
    #similar to __init__
    def add_run(self, **kwargs):
        run = self.run
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
        self.step_n_jumps.append(kwargs.get('step_n_jumps', int(self.n_steps[run] / 1000)))
        self.step_min_val.append(kwargs.get('step_min_val', 0.))
        self.step_max_val.append(kwargs.get('step_max_val', 100.))

        self.values_p.append(kwargs.get('values_p', self.values_p[self.run]))
        self.probabilities_p.append(kwargs.get('probabilities_p', self.probabilities_p[run]))
        self.s.append(kwargs.get('s', h.random_cont(self.n_heads[run], -self.d[self.run]/2, self.d[self.run]/2)))
        self.p.append(kwargs.get('p', h.random_discrete(self.n_heads[run], probabilities=self.probabilities_p[self.run+1])))
        self.t0.append(kwargs.get('t0', self.t0[self.run]))

        self.run += 1
        #OBSOLETE if self.v_Coeff[self.run] != 0 and self.mode[self.run] == 'fControl': raise ValueError('Wrong mode or velocity parameters chosen')
        run = self.run

        if self.mode[self.run] in ['vControl']: self.t.append(np.array([i*self.d_t[self.run] + self.t0[self.run] for i in range(self.n_steps[self.run])]))
        elif self.mode[self.run] in ['fControl']: self.t.append(np.zeros(self.n_steps[run]))
        self.sum_P.append(np.array([]))
        self.n_neighbours.append(h.find_neighbours(self.th[self.run], self.d[self.run], self.bta[self.run], self.k[self.run], self.k_on[run]))


        #append new slots for storage of desired Variables out of s,p,f,sum_F
        if self.store.get('s'): self.S.append(np.zeros((self.n_steps[run], self.n_heads[run])))
        if self.store.get('pos') or self.mode[self.run] in ['fControl', 'springControl']: self.Pos.append(np.zeros((self.n_steps[self.run],1)))
        if self.store.get('pos_pull') or self.mode[self.run] in ['springControl']: self.Pos_pull.append(np.zeros((self.n_steps[self.run],1)))
        if self.store.get('p'): self.P.append(np.zeros((self.n_steps[run], self.n_heads[run])))
        if self.store.get('f'): self.F.append(np.zeros((self.n_steps[run], self.n_heads[run])))
        if self.store.get('sum_f'): self.sum_F.append(np.zeros((self.n_steps[run],1)))


        self.f_Sum_axis.append([])
        self.v_axis.append([])

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
    def update_vC(self, s, p, r01, run, step, bta, k):
        #update position according to option
        if self.option[run] == 'poly':
            displ = np.polynomial.polynomial.polyval(step, self.v_Coeff[run])
        elif self.option[run] == 'step':
            displ = h.step_fcn(self.n_steps[run], self.step_n_jumps[run], self.step_min_val[run], self.step_max_val[run], step)
        elif self.option[run] == 'const':
            displ = self.v_Coeff[run][0]
        else:
            raise ValueError("error! option: {}".format(self.option[run]))
        s = h.translate_def(s, displ * self.d_t[run])

        #case: unbound
        if h.unitize(p) == 0:
            K_plus = h.k_plus(h.s_row(h.wrapping(s, self.d[run]), self.n_neighbours[run], self.d[run]), p, self.d[run], self.bta[run], self.k[run], self.k_on[run])
            K_sum = sum(K_plus)
            if r01 <= sum(K_plus) * self.d_t[run]:
                p = h.det_p(h.p_row(self.n_neighbours[run]), K_plus, K_sum)
                s = h.wrapping(s, self.d[run])
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
    
    def update_sC(self, s, p, d, K_plus, K_sum, rn):
        
    #MAIN UPDATE METHODS#END############################################

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
        ax.set_ylabel('Total load')

        for r in run:
            if len(self.args_passed[r]) == 1:
                lab = list(self.args_passed[r])[0]
                val = self.args_passed[r][lab]
                ax.plot(self.t[r], self.sum_F[r] / float(self.n_heads[r]), linewidth=1.0, linestyle="-", label='{}={}'.format(lab, val))
            else: ax.plot(self.t[r], self.sum_F[r] / float(self.n_heads[r]), linewidth=1.0, linestyle="-")
        if leg: ax.legend()
        plt.savefig('Sum_f.png', dpi=200)

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

    #run actual simulation
    def start_run(self, run):
        #don't touch the object variable
        s = self.s[run]
        p = self.p[run]
        sum_F = 0.
        pos = 0.
        t = 0.

        #~~~~~~~~~~~~~~~~Here, the seed of the random generator could be implemented~~~~~~~~~~~~~~~~~~~~~~
        #random numbers for udating (for 'fControl', we need two sets)
        if self.mode[run] in ['vControl']:
            print('rand creates')
            rand01 = random.random_sample((self.n_steps[run], self.n_heads[run]))
        if self.mode[run] in ['fControl']:
            rand011 = random.random_sample(self.n_steps[run])
            rand012 = random.random_sample(self.n_steps[run])
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
            rand011 = random.random_sample(self.n_steps[run])
            rand012 = random.random_sample(self.n_steps[run])
            rand013 = random.random_sample(self.n_steps[run])
            
            #stretching to desired force
            f = h.forceV(s, p, self.d[run])

            sum_F = sum(f)
            #~ print sum_F
            n_att = sum(h.unitizeV(p))
            #~ print n_att
            displ = - sum_F / (n_att + self.k_pull[run])
            #~ print displ
            s = h.translateV(s, displ)
            pos += displ

        print('having generated rands')
        print(pos)
        #loop over all iteration steps
        for i in tqdm.tqdm(range(self.n_steps[run])):
            #show, how many iterations have already elapsed
            #~ if float(i) / self.n_heads[run] in [.1, .2, .5, 1.]:
                #~ print str(float(i) / self.n_heads[run] *100) + '% of steps elapsed (run {})'.format(run)

            #calculate force, f contains all forces of each head, sum_F is total load on filament
            f = h.forceV(s, p, self.d[run])
            sum_F = sum(f)

            #store Variables s, p
            if self.store.get('s'): self.S[run][i] = s
            if self.store.get('pos'): self.Pos[run][i] = pos
            if self.store.get('p'): self.P[run][i] = p
            if self.store.get('f'): self.F[run][i] = f
            if self.store.get('sum_f'): self.sum_F[run][i] = sum_F


            #update s and p according to mode
            #self.s and self.p remain unchanged!
            if self.mode[run] == 'vControl':
                s, p = self.updateV_vC(self, s, p, rand01[i], run, i, self.bta[run], self.k[run])
                if self.option[run] == 'poly':
                    displ = np.polynomial.polynomial.polyval(i, self.v_Coeff[run])
                elif self.option[run] == 'step':
                    displ = h.step_fcn(self.n_steps[run], self.step_n_jumps[run], self.step_min_val[run], self.step_max_val[run], i)
                pos += displ
            elif self.mode[run] == 'fControl':
                #here, probabilities for attaching, detaching and the corresponding wating time tau are calculated
#                s_mat_w = h.s_matrix(h.wrapping(s, self.d[run]), self.n_neighbours[run], self.d[run]) 
#                k_plus_mat = h.k_plusV(s_mat_w, p.reshape(self.n_heads[run],1), self.d[run], self.bta[run], self.k[run], self.k_on[run])
#                k_plus_sum = k_plus_mat.sum(axis=1)
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
                
                f = h.forceV(s, p, self.d[run])

                f_delta = sum(f) - self.loadF[run]


                #calculate number of attached heads
                n_att = sum(h.unitizeV(p))
                #calculate force <-> displacement of filament
                if n_att == 0:
                    print(i)
                    print("connection broke!")
                    #~ for j in [l+i for l in range(self.n_steps[run] - i)]:

                    break

                displ = -f_delta / float(n_att)

                #updating filament position and s positions
                pos += displ
                s = h.translateV(s, displ)

                #updating the elapsed time
                #~ print tau_min
                t += tau
                self.t[run][i] = t
                
            elif self.mode[run] == 'springControl':
                

        #preparing the time vector to be added in first column of each file
        t_w = self.t[run][np.newaxis]
        t_w = t_w.T
        #~ print t_w
        #~ print self.Pos[run]
        #writing the results to textfiles
        if self.writeText:
            if self.store.get('s'):
                print('Store S')
                S_out = np.concatenate((t_w, self.S[run]), axis=1)
                np.savetxt('S_{}.txt'.format(run), S_out, header='time, Position of each head')
            if self.store.get('pos'):
                print('Store Pos')
                Pos_out = np.concatenate((t_w, self.Pos[run]), axis=1)
                np.savetxt('Pos_{}.txt'.format(run), Pos_out, header='time, Position of filament')
            if self.store.get('p'):
                print('Store P')
                P_out = np.concatenate((t_w, self.P[run]), axis=1)
                np.savetxt('P_{}.txt'.format(run), P_out, header='time, binding state of each head')
            if self.store.get('f'):
                print('Store F')
                F_out = np.concatenate((t_w, self.F[run]), axis=1)
                np.savetxt('F_{}.txt'.format(run), F_out, header='time, force applied by each head')
            if self.store.get('sum_f'):
                sum_F_out = np.concatenate((t_w, self.sum_F[run]), axis=1)
                np.savetxt('sum_F_{}.txt'.format(run), sum_F_out, header='time, force applied by each head')

        return s, p

    def create_threads(self, runs):
        self.jobs = [th.Thread(target=self.start_run, args=(i,)) for i in runs]
    def start_threads(self):
        for j in self.jobs:
            j.start()
        for j in self.jobs:
            j.join()

	#get the asymptotic value for each step of v_step
    def average_norm_force(self, run, f_Sum, v_step_list, equilib_wait_frac=0.2):

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
        #~ print(self.f_Sum_axis_added)
        #~ print(self.v_axis_added)

#initialization

#################################################################################################|

mode = 'fControl' #choose from ['vControl', 'fControl', ]
option = 'poly' #for fControl choose from xy			
				#for vControl choose from 													
				#						-> poly: specify coefficients						
				#						-> step: specify n_elem, n_jumps, min_val, max_val		
name = 'matematica_crosscheck_f=1'													
																								
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
loadF = [1. for i in range(3)]																
																								
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