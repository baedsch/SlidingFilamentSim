import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
from numpy import random
#~ import random
import math
import os
import time as tme
import time
import threading as th
import multiprocessing as mp


# in order to be able to use map with KeyWordArgs
from functools import partial
import helper as h




class simulation:
	#initialization:
	#kwargs: mode, seed, d, s, p, d_t, t0, k, bta, values_p, probabilities_p, vCoefficients (function of v polynomial)
	#checks, which variables should be stored, by default, everything is stored
	def __init__(self, mode, n_steps=int(1e4), n_heads=int(1e3), name="default", **kwargs):
		#~ multiprocessing.Process.__init__(self)
		#these are the parameters
		#self.variables GET CHANGED WHITHIN FCN SCOPES!! no globaliyation needed (self. is like global)
		#most values are stored np.arrays nested in lists: for each run, append new array to list
		#within the arrays, axis 0 contains the time, axis 1 stands for the head index
		self.name = str(name)
		self.run = 0
		self.n_steps = [n_steps]
		self.n_heads = [n_heads]
		if mode in ['vControl', 'fControl', ]: self.mode = [mode]
		else: raise ValueError('Wrong mode chosen')
		self.option = [kwargs.get('option', '')]
		if not kwargs.get('option', '') in ['poly', 'step', '']: raise ValueError('Wrong option chosen')
		self.args_passed = [kwargs]
		
		self.d = [kwargs.get('d', 1.)]
		self.d_t = [kwargs.get('d_t', 1e-3)]
		self.seed = [kwargs.get('seed', False)]
		
		self.loadF = [kwargs.get('loadF', 200.)]
		self.bta = [kwargs.get('bta', 2.)]
		self.k = [kwargs.get('k', 2.)]
		self.th = [kwargs.get('th', 0.0001)]
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
		
		self.n_neighbours = [h.find_neighbours(self.th[self.run], self.d[self.run], self.bta[self.run], self.k[self.run])]
		
		#declare storage of desired Variables out of s,p,f,sum_F
		self.store = {}
		self.store['s'] = kwargs.get('s_store',True)
		if self.store.get('s'): self.S = [np.zeros((n_steps, n_heads))]
		self.store['pos'] = kwargs.get('pos_store',True)
		if self.store.get('pos') or self.mode[self.run] == 'fControl': self.Pos = [np.zeros((n_steps,1))]
		self.store['p'] = kwargs.get('p_store',True)
		if self.store.get('p'): self.P = [np.zeros((n_steps, n_heads))]
		self.store['f'] = kwargs.get('f_store',True)
		if self.store.get('f'): self.F = [np.zeros((n_steps, n_heads))]
		self.store['sum_f'] = kwargs.get('sum_f_store',True)
		if self.store.get('sum_f'): self.sum_F = [np.zeros((n_steps,1))]
		
		#write the stored variables into textfiles?
		self.writeText = kwargs.get('writeText', False)
		#make directory
		directory = tme.strftime("%Y_%m_%d__%H:%M:%S_", tme.gmtime()) + self.name
		try:
			os.mkdir('./res/' + directory)
			os.chdir('res/' + directory)
		except:
			print "Directory already existing, please choose unique simulation name! Otherwise data might get overwritten"
		
	#similar to __init__
	def add_run(self, **kwargs):
		run = self.run
		self.n_steps.append(kwargs.get('n_steps', self.n_steps[self.run]))
		self.n_heads.append(kwargs.get('n_heads', self.n_heads[self.run]))
		mode = kwargs.get('mode', self.mode[self.run])
		if mode in ['vControl', 'fControl', ]: self.mode.append(kwargs.get('mode', self.mode[self.run]))
		else: raise ValueError('Wrong mode chosen')
		self.option.append(kwargs.get('option', ''))
		if not kwargs.get('option', '') in ['poly', 'step', '']: raise ValueError('Wrong option chosen')
		self.args_passed.append(kwargs)
		
		self.d.append(kwargs.get('d', self.d[self.run]))
		self.d_t.append(kwargs.get('d_t', self.d_t[self.run]))
		self.seed.append(kwargs.get('seed', self.seed[self.run]))
		
		self.loadF.append(kwargs.get('loadF', self.loadF[self.run]))
		self.bta.append(kwargs.get('bta', self.bta[self.run]))
		self.k.append(kwargs.get('k', self.k[self.run]))
		self.th.append(kwargs.get('th', self.th[self.run]))
		self.v_Coeff.append(kwargs.get('v_Coeff', self.v_Coeff[self.run]))
		self.step_n_jumps.append(kwargs.get('step_n_jumps', int(n_steps / 1000)))
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
		self.n_neighbours.append(h.find_neighbours(self.th[self.run], self.d[self.run], self.bta[self.run], self.k[self.run]))
		
		
		#append new slots for storage of desired Variables out of s,p,f,sum_F
		if self.store.get('s'): self.S.append(np.zeros((self.n_steps[run], self.n_heads[run])))
		if self.store.get('pos') or self.mode[self.run] == 'fControl': self.Pos.append(np.zeros((self.n_steps[self.run],1)))
		if self.store.get('p'): self.P.append(np.zeros((self.n_steps[run], self.n_heads[run])))
		if self.store.get('f'): self.F.append(np.zeros((self.n_steps[run], self.n_heads[run])))
		if self.store.get('sum_f'): self.sum_F.append(np.zeros((self.n_steps[run],1)))
		
	# initialize a randomly chosen s (resp. p) as a startvector
	def init_s_rand(self):
		self.s[self.run] = random_cont(self.n_heads, -self.d/2, self.d/2)
		return self.s[self.run]
	def init_p_rand(self):
		self.p[self.run] = random_discrete(self.n_heads, self.values_p[self.run], self.probabilities_p[self.run])
		return self.p[self.run]
		
	#################################################MAIN UPDATE METHODS
	#main update method for mode==vControl
	def update_vC(self, s, p, r01, run, step, bta, k):
		#update position according to option
		if self.option[run] == 'poly':
			displ = np.polynomial.polynomial.polyval(step, self.v_Coeff[run])
		elif self.option[run] == 'step':
			displ = h.step_fcn(self.n_steps[run], self.step_n_jumps[run], self.step_min_val[run], self.step_max_val[run], step)
		else: 
			print "error!!!!!!!!!!!!!!!!!!!!111, option:"
			print self.option[run]
		s = h.translate_def(s, displ)
		
		#case: unbound
		if h.unitize(p) == 0:
			K_plus = h.k_plus(h.s_row(h.wrapping(s, self.d[run]), self.n_neighbours[run], self.d[run]), p, self.d[run], self.bta[run], self.k[run])
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
	def update_fC(self, s, p, K_plus, K_sum, run):
		#case: unbound
		if h.unitize(p) == 0:
			p = h.det_p(h.p_row(self.n_neighbours[run]), K_plus, K_sum)
			s = h.wrapping(s, self.d[run])
		#case: bound
		else:
			p = 0
			
		return s, p
	#vectorize it in order to enhance performance
	updateV_fC = np.vectorize(update_fC)
	
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
			print 'input into list'
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		
		ax.set_title('Time evolution of bound heads ({} run(s)), mode: {}'.format(len(run), self.mode[run[0]]))
		ax.set_xlabel('time [s]')
		ax.set_ylabel('bound heads')
		
		for r in run:
			if len(self.args_passed[r]) == 1:
				lab = self.args_passed[r].keys()[0]
				val = self.args_passed[r][lab]
				ax.plot(self.t[r], self.sum_P[r], linewidth=1.0, linestyle="-", label='{}={}'.format(lab, val))
			else: ax.plot(self.t[r], self.sum_P[r], linewidth=1.0, linestyle="-")
		if leg: legend = ax.legend()
		plt.savefig('Sum_p.png', dpi=200)
	
	def plot_f(self, run, leg=False):
		if isinstance(run, int):
			run = [run]
			print 'input into list'
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		
		ax.set_title('Total load on filament ({} run(s)), mode: {}'.format(len(run), self.mode[run[0]]))
		ax.set_xlabel('time [s]')
		ax.set_ylabel('Total load')
		
		for r in run:
			if len(self.args_passed[r]) == 1:
				lab = self.args_passed[r].keys()[0]
				val = self.args_passed[r][lab]
				ax.plot(self.t[r], self.sum_F[r], linewidth=1.0, linestyle="-", label='{}={}'.format(lab, val))
			else: ax.plot(self.t[r], self.sum_F[r], linewidth=1.0, linestyle="-")
		if leg: legend = ax.legend()
		plt.savefig('Sum_f.png', dpi=200)
	
	def plot_pos(self, run, leg=False):
		if isinstance(run, int):
			run = [run]
			print 'input into list'
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		
		ax.set_title('Position filament ({} run(s)), mode: {}'.format(len(run), self.mode[run[0]]))
		ax.set_xlabel('time [s]')
		ax.set_ylabel('position')
		
		for r in run:
			if len(self.args_passed[r]) == 1:	
				lab = self.args_passed[r].keys()[0]
				val = self.args_passed[r][lab]
				ax.plot(self.t[r], self.Pos[r], linewidth=1.0, linestyle="-", label='{}={}'.format(lab, val))
			else: ax.plot(self.t[r], self.Pos[r], linewidth=1.0, linestyle="-")
		if leg: legend = ax.legend()
		plt.savefig('Pos.png', dpi=200)
	
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
			print 'rand creates'
			rand01 = random.random_sample((self.n_steps[run], self.n_heads[run]))
		if self.mode[run] in ['fControl']: 
			rand01_plus = random.random_sample((self.n_steps[run], self.n_heads[run]))
			rand01_min = random.random_sample((self.n_steps[run], self.n_heads[run]))
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
			
		print 'having generated rands'
		print pos
		#loop over all iteration steps
		for i in range(self.n_steps[run]):
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
				s_mat = h.s_matrix(h.wrapping(s, self.d[run]), self.n_neighbours[run], self.d[run])
				k_plus_mat = h.k_plusV(s_mat, p.reshape(self.n_heads[run],1), self.d[run], self.bta[run], self.k[run])
				k_plus_sum = k_plus_mat.sum(axis=1)
				k_min = h.k_minV(s, p, self.bta[run], self.k[run])
				tau_p = np.array(h.tauV(rand01_plus[i], k_plus_sum))
				
				tau_m = np.array(h.tauV(rand01_min[i], k_min))
				tau = [plus+minus for plus,minus in zip(tau_p, tau_m)]
				#find the result for the waiting time
				tau_min = np.min(tau)
				#find head with lowest waiting time, if more heads have the same 
				min_index = np.argwhere(tau == tau_min)[0]
				
				if len(min_index) > 1: index = random.randint(len(min_index))
				else: index = 0
				min_index = min_index[index]
				
				#update this head, calculate the force difference
				s_i, p_i = s[min_index], p[min_index]
				f_i = h.force(s_i, p_i, self.d[run])
					#in case of binding, s gets wrapped
				s_upd, p_upd = self.update_fC(s_i, p_i, k_plus_mat[min_index], k_plus_sum[min_index], run)
				s[min_index] = s_upd
				p[min_index] = p_upd
				#~ print (p_i, p_upd, min_index, tau_min)
				f_upd = h.force(s_upd, p_upd, self.d[run])
				f_delta = f_upd - f_i
				
				#calculate number of attached heads
				n_att = sum(h.unitizeV(p))
				#calculate force <-> displacement of filament
				if n_att == 0: 
					print i
					print "connection broke!"
					#~ for j in [l+i for l in range(self.n_steps[run] - i)]:
						
					break
					
				displ = -f_delta / float(n_att)
				
				#updating filament position and s positions
				pos += displ
				s = h.translateV(s, displ)
				
				#updating the elapsed time
				#~ print tau_min
				t += tau_min
				self.t[run][i] = t
				
		#preparing the time vector to be added in first column of each file
		t_w = self.t[run][np.newaxis]
		t_w = t_w.T
		#~ print t_w
		#~ print self.Pos[run]
		#writing the results to textfiles
		if self.writeText:
			if self.store.get('s'):
				print 'Store S'
				S_out = np.concatenate((t_w, self.S[run]), axis=1)
				np.savetxt('S_{}.txt'.format(run), S_out, header='time, Position of each head')
			if self.store.get('pos'):
				print 'Store Pos'
				Pos_out = np.concatenate((t_w, self.Pos[run]), axis=1)
				np.savetxt('Pos_{}.txt'.format(run), Pos_out, header='time, Position of filament')
			if self.store.get('p'):
				print 'Store P'
				P_out = np.concatenate((t_w, self.P[run]), axis=1)
				np.savetxt('P_{}.txt'.format(run), P_out, header='time, binding state of each head')
			if self.store.get('f'):
				print 'Store F'
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

	#~ def create_threadsMP(self, runs):
		#~ self.jobs = [mp.Process(target=self.start_run, args=(i,)) for i in runs]
	#~ def start_threadsMP(self):
		#~ for j in self.jobs:
			#~ j.start()
		#~ for j in self.jobs:
			#~ j.join()

#~ class myth(threading.Thread):
	#~ def run(self):
		#~ global		
#~ $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$			
#~ from multiprocessing import Process, Queue

#~ def f(q):
    #~ q.put([42, None, 'hello'])

#~ if __name__ == '__main__':
    #~ q = Queue()
    #~ p = Process(target=f, args=(q,))
    #~ p.start()
    #~ print q.get()    # prints "[42, None, 'hello']"
    #~ p.join()
#~ $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$https://docs.python.org/2/library/multiprocessing.html#the-process-class
#~ queues + put is the spirit (probably)
