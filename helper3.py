import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
from numpy import random
#~ import random
import os
import time as tme
import time


#wrapping fcn
def wrapping(x, d):
    return (x - d / 2) % d - d / 2
wrappingV = np.vectorize(wrapping)

#convert vector elementwise in 0. if eq. 0 or 1. else
def unitize(array):
    if type(array).__module__ == np.__name__:
        arrayUni = array != 0
        return arrayUni.astype(float)
    if type(array).__name__ in ['int', 'float']:
        return float(array != 0)
    else:
        raise AssertionError('Has to be numpy array, int or float')
#vectorize it in order to enhance performance
unitizeV = np.vectorize(unitize)

#force of each head (also the delta_s between the head and the binding site), to be used with map fcn
def force(s, p, d):
    #case head not attached
    if p == 0:
        return 0.
    #case head position in closest interval or s >= d/2
    if p > 0:
        return s + (p - 1) * d
    #case head position  s <= -d/2
    if p < 0:
        return s + p * d
#vectorize it in order to enhance performance
forceV = np.vectorize(force)

#detaching rate(s,p,beta,k) to be used with map fcn
def k_min(s, p, bta, k):
    pU = unitize(p)
    sSq = s**2

    if abs(s) < (1 + k) / 2:
        return pU * np.exp(-(bta * sSq) / (1 + k)) * np.cosh(bta * s)

    else:
        return 100000000000000
k_minV = np.vectorize(k_min)

#attaching rate (+V) to be used with map fcn, either feed already wrapped s in or set w=True!
def k_plus(s, p, d, bta, k, k_on,  w=False):
    if w: s = wrapping(s, d)
    pU = unitize(p)
    b = np.sqrt(bta)
    c = .5
    #solution of gaussian integral
    gaussianInt = .5 * k_on * ( - special.erf(b * (s - c)) + special.erf(b * (s + c)))

    return (1 - pU) * gaussianInt
k_plusV = np.vectorize(k_plus)

#attaching rate sum, either feed already wrapped s in or set w=True!
def k_plus_sum(s, p, d, bta, k, k_on, n_neighbours,  w=False):
    if w: s = wrapping(s, d)
    pU = unitize(p)
    b = np.sqrt(bta)
    c = .5
    res = 0
    for i in [- n_neighbours + z for z in range(2 * n_neighbours + 1)]:
        #solution of gaussian integral
        res += .5 * k_on * ( - special.erf(b * (s + (i * d) - c)) + special.erf(b * (s + + (i * d) + c)))

    return (1 - pU) * res
k_plus_sumV = np.vectorize(k_plus_sum)

#returns the waiting time until reaction occurs
def tau(r, k):
    if k == 0: return 0.
    else: return - np.log(r) / k
tauV = np.vectorize(tau)

#translate
def translate_def(s, dist):
    return s + dist

translateV = np.vectorize(translate_def)

#helper fcn
def gaussian(x, a, b):
    return b * np.exp(-a * x**2)

#returns image of stepfunction determined by n_elem, n_jumps, min_val, max_val
def step_list(n_elem, n_jumps, min_val, max_val):
    c = [int(i / (n_elem / n_jumps)) for i in range(n_elem)]
    a = [i * (max_val - min_val) / float(n_jumps - 1) + min_val for i in c]
    return a

#returns element of stepfunction determined by n_elem, n_jumps, min_val, max_val
def step_fcn(n_elem, n_jumps, min_val, max_val, index):
    c = int(index / (n_elem / n_jumps))
    a = c * (max_val - min_val) / float(n_jumps - 1) + min_val
    return a

def detach_plot_fcn(s, bta, k):
    result = []
    for s_iter in s:

        if np.absolute(s_iter) < (1 + k) / 2.:
            result.append( np.exp(-(bta * s_iter**2) / (1. + k)) * np.cosh(bta * s_iter))

        else:
            result.append(  1.)

    return result

def plot_detach():
    bta = 0.6
    k = 10.
    X = np.linspace(-10., 10., 256)
    D = detach_plot_fcn(X, bta, k)
    plt.figure(figsize=(8,6), dpi=80)
    plt.subplot(111)
    plt.plot(X, D, color="blue", linewidth=1.0, linestyle="-")

    plt.show()

def plot_attach():
    bta = 2.
    k = 10
    k_on  = 10.
    X = np.linspace(-1., 1., 2560)
    D = k_plus_sum(X, 0, 2, bta, k, k_on,5,)
    plt.figure(figsize=(8,6), dpi=80)
    plt.subplot(111)
    plt.plot(X, D, color="blue", linewidth=1.0, linestyle="-")

    plt.show()
plot_attach()
#find the number of relevant neighbours up to arbitrary theshold
def find_neighbours(th, d, bta, k, k_on):
    p = 0

    n = 0
    while k_plus(d * (n + .5), p, d, bta, k, k_on, w=False) / k_plus(0, p, d, bta, k, k_on, w=False) > th:
        n += 1
    return n

#returns 2dArray with rows as following [s - (n_neighbours * d), ..., s, ..., s + (n_neighbours * d)]
def s_matrix(s, n_neighbours, d):
    n_heads = len(s)

    #2n + 1 for dim of matrix
    neighbours_ext = 2 * n_neighbours + 1
    #dimension of s matrix
    dim = (n_heads, neighbours_ext)
    #returns 2 dim array acc to dim filled with float 0.
    s_mat = np.zeros(dim, dtype=float)

    #s_matrix with given s-values or randomly assign
    if isinstance(s, int):
        random.seed(1234789)
        s = [d * random.random() - d / 2 for i in range(n_heads)]
    else:
        n_heads = len(s)

    #assigning the matrix
    for r,i in zip(s,range(n_heads)):
        for j in range(neighbours_ext):
            s_mat[i][j] = r - n_neighbours * d + j * d

    return s_mat

#returns row as following [s - (n_neighbours * d), ..., s, ..., s + (n_neighbours * d)]
def s_row(s, n_neighbours, d):

    #2n + 1 for dim of matrix
    neighbours_ext = 2 * n_neighbours + 1
    #dimension of s matrix

    #returns 2 dim array acc to dim filled with float 0.
    s_mat = np.zeros(neighbours_ext, dtype=float)

    #assigning the row
    for j in range(neighbours_ext):
        s_mat[j] = s - n_neighbours * d + j * d

    return s_mat

def p_row(n_neighbours):

    #2n + 1 for dim of matrix
    neighbours_ext = 2 * n_neighbours + 1
    #dimension of s matrix

    #returns array filled with int 0
    p_mat = np.zeros(neighbours_ext, dtype=int)

    #assigning the row
    i = 1
    for j in [n_neighbours + n for n in range(n_neighbours + 1)]:
        p_mat[j] = i
        i += 1
    i = -1
    for j in [n_neighbours - 1 - n for n in range(n_neighbours)]:
        p_mat[j] = i
        i -= 1
    return p_mat

#there is a function which does this np.random.choise(values, p=probabilities)
#could be improved by letting the fcn eat pregenerated rand
def det_p(values, probabilities, norm):
    probabilities = [p / norm for p in probabilities]
    if not len(values) == len(probabilities):
        raise AssertionError('arg1 and arg2 must have the same length')

    #returns e.g. [0., 0.3, 1.0], therefore takes element and adds the following one
    acc = np.add.accumulate(probabilities)
    bins = np.append(0., acc)
    r = random.random()
    for i in range(len(values)):
        if bins[i] <= r < bins[i+1]: return values[i]
    raise AssertionError('something in the det_p method went WRONG!!!')


#returns 2dArray with attachment rates as entries acc to s_matrix values (see above) OBSOLETE
def k_plus_matrix(s_mat, p, d, bta, k, k_on):
    return     k_plus(s_mat, p, d, bta, k, k_on)
k_plus_matrixV = np.vectorize(k_plus_matrix)

#there is a function which does this np.random.choise(values, p=probabilities)
def random_discrete(length, values=[0., 1.], probabilities=[0.3, 0.7]):
    if not np.sum(probabilities) == 1.:
        raise AssertionError('Probs should add to 1; input dicts should have the form length, values=[0., 1.], probabilities=[0.3, 0.7]')
    if not len(values) == len(probabilities):
        raise AssertionError('arg1 and arg2 must have the same length')
    res = np.array([])

    #returns e.g. [0., 0.3, 1.0], therefore takes element and adds the following one
    acc = np.add.accumulate(probabilities)
    bins = np.append(0., acc)
    for n in range(length):
        r = random.random()
        for i in range(len(values)):
            if bins[i] <= r < bins[i+1]: res = np.append(res, values[i])
    return res

#there might be also a lib fcn for this
def random_cont(length, lower, upper):
    res = np.array([])
    interval = float(upper - lower)
    for n in range(length):
        res = np.append(res, interval * random.random() - interval / 2.)
    return res
