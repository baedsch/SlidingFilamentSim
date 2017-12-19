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
        raise ValueError('Has to be numpy array, int or float')

#vectorize it in order to enhance performance
unitizeV = np.vectorize(unitize)

#force of each head (also the delta_s between the head and the binding site), to be used with map fcn
def force(s, p, d):
    #s ~y for attached heads, see update functions
    p = p != 0
    #elementwise
    return s * p

#vectorize it in order to enhance performance
forceV = np.vectorize(force)

#detaching rate(s,p,beta,k) to be used with map fcn
def k_min(s, p, bta, k):
    #s ~y for attached heads, see update functions
    pU = unitize(p)
    sSq = s**2

    if abs(s) < (1 + k) / 2:
        return pU * np.exp(-(bta * sSq) / (1 + k)) * np.cosh(bta * s)

    else:
        return 1e50
k_minV = np.vectorize(k_min)

#integral of k_min
def int_k_min(s1, s2, p, bta, k, v):
    root = np.sqrt(bta / (1 + k))

    #                                              -V- this v added brcause of variable change
    c = np.sqrt(np.pi * (k + 1) / bta) / (4 * bta * v)
    intt = lambda s: c * (np.exp(bta * (1 + k) / 4) * (special.erf(root * (1 + k + 2*s) / 2) - special.erf(root * (1 + k - 2*s) / 2)))

    return abs(intt(s2) - intt(s1))

#attaching rate (+V) to be used with map fcn, either feed already wrapped s in or set w=True!
def k_plus(s, p, d, bta, k, k_on,  w=True):
    if w: s = wrapping(s, d)
    pU = unitize(p)
    b = np.sqrt(bta)
    c = .5
    #solution of gaussian integral
    gaussianInt = .5 * k_on * ( - special.erf(b * (s - c)) + special.erf(b * (s + c)))

    return (1 - pU) * gaussianInt
k_plusV = np.vectorize(k_plus)

#the lambda expression
def int_k_plus_sum(s1, s2, d, bta, k, k_on, n_neighbours, v, w=False):
    #this is andatory because direction of integration does not matter MEESES UP!
#    if s2 < s1: s1, s2 = s2, s1
#    print('s1, s2.................................................................',s1,'...',s2)
    s1t = s1 + d/2
    s1tw = wrapping(s1t, d)
    s2t = s2 + d/2
    s2tw = wrapping(s2t, d)

    b = np.sqrt(bta)
    spi = np.sqrt(np.pi)
    c = 0.5

    integralt = lambda s, i : c * k_on / (b * spi * v) * (np.exp(-bta * ((c + i*d + (s-d/2))**2)) - np.exp(-bta * ((-c + i*d + (s-d/2))**2))) + 1/v * c * k_on * ((c + i*d + (s-d/2)) * special.erf(b * (c + i*d + (s-d/2))) - (-c + i*d + (s-d/2)) * special.erf(b * (-c + i*d + (s-d/2))))

    int_full_period = lambda i: integralt(0, i) - integralt(d, i)

    if s1 == s2: return 0.
    if s2tw >= s1tw and s2 >= s1:
        res = 0
        for n in [i - n_neighbours for i in range(2 * n_neighbours + 1)]:

            #AAAACHtung: n oder 0 in lambda???
            res += integralt(s2tw, n) - integralt(s1tw, n) + (abs((s2t - s1t)) // d) * int_full_period(n) #// is floor division in python
#        print('res.........................................................',res)
        return res
    elif s2tw < s1tw and s2 >= s1:
        res = 0
        for n in [i - n_neighbours for i in range(2 * n_neighbours + 1)]:
           res += integralt(s2tw, n) - integralt(s1tw, n) + (abs((s2t - s1t)) // d + 1) * int_full_period(n) #// is floor division in python
#        print('res.........................................................',res)
        return res
    elif s2tw >= s1tw and s2 < s1:
        res = 0
        for n in [i - n_neighbours for i in range(2 * n_neighbours + 1)]:

            #AAAACHtung: n oder 0 in lambda???
            res += integralt(s2tw, n) - integralt(s1tw, n) + (abs((s2t - s1t)) // d + 1) * int_full_period(n) #// is floor division in python
#        print('res.........................................................',res)
        return res
    elif s2tw < s1tw and s2 < s1:
        res = 0
        for n in [i - n_neighbours for i in range(2 * n_neighbours + 1)]:
           res += integralt(s2tw, n) - integralt(s1tw, n) + (abs((s2t - s1t)) // d) * int_full_period(n) #// is floor division in python
#        print('res.........................................................',res)
        return res
    else: raise ValueError()

#attaching rate sum, either feed already wrapped s in or set w=True!
def k_plus_sum(s, p, d, bta, k, k_on, n_neighbours,  w=True):
    if w: s = wrapping(s, d)
    pU = unitize(p)
    b = np.sqrt(bta)
    c = .5
    res = 0
    for i in [- n_neighbours + z for z in range(2 * n_neighbours + 1)]:
        #solution of gaussian integral
        res += .5 * k_on * ( - special.erf(b * (s + (i * d) - c)) + special.erf(b * (s + (i * d) + c)))

    return (1 - pU) * res
k_plus_sumV = np.vectorize(k_plus_sum)

#get minimum and maximum value of k_plus_sum within subintervals
def min_max_plus(n, d, bta, k, k_on, n_neighbours):
    res = {}
    varialbe_values = [- d / 2 + z * d / (2 * n) for z in range(2 * n + 1)]
    f_values = [k_plus_sum(i, 0, d, bta, k, k_on, n_neighbours) for i in varialbe_values]
    for i in range(len(varialbe_values) - 1):
        res[varialbe_values[i]] = sorted([f_values[i], f_values[i + 1]])

    return res

#get minimum and maximum value of k_min within subintervals
def min_max_min(n, bta, k):
    res = {}
    d = (1 + k)
    varialbe_values = [- d / 2 + z * d / (2 * n) for z in range(2 * n + 1)]
    f_values = [k_min(i, 1, bta, k) for i in varialbe_values]
    for i in range(len(varialbe_values) - 1):
        res[varialbe_values[i]] = sorted([f_values[i], f_values[i + 1]])

    return res

def get_max_k_min(bta, k):
    return np.exp(-(bta * ((1 + k) / 2)**2 ) / (1 + k)) * np.cosh(bta * (1 + k) / 2)

def get_max_k_plus_sum(d, bta, k, k_on, n_neighbours):
    b = np.sqrt(bta)
    c = .5
    res = 0
    for i in [- n_neighbours + z for z in range(2 * n_neighbours + 1)]:
        #solution of gaussian integral
        res += .5 * k_on * ( - special.erf(b * (0 + (i * d) - c)) + special.erf(b * (0 + + (i * d) + c)))

    return res
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
            result.append(  2.7)

    return result

def plot_detach():
    bta = 2
    k = 10.
    X = np.linspace(-10., 10., 256)
    D = detach_plot_fcn(X, bta, k)
    plt.figure(figsize=(8,6), dpi=80)
    plt.subplot(111)
    plt.plot(X, D, color="blue", linewidth=1.0, linestyle="-")

    plt.show()
#plot_detach()


def plot_attach_diff():
    bta = 2.
    k = 10
    k_on  = 10.
    X = np.linspace(-3.,3., 2560)
    D = k_plus_sum(X, 0, 0.2, bta, k, k_on,10,w=True) - k_plus_sum(X, 0, 0.2, bta, k, k_on,9,w=True)
    plt.figure(figsize=(8,6), dpi=80)
    plt.subplot(111)
    plt.plot(X, D, color="blue", linewidth=1.0, linestyle="-")

    plt.show()
#plot_attach_diff()

def plot_attach_one():
    bta = 2.
    k = 10
    k_on  = 10.
    X = np.linspace(-3.,3., 2560)
    D = k_plus_sum(X, 0, 0.2, bta, k, k_on,2,w=False)
    plt.figure(figsize=(8,6), dpi=80)
    plt.subplot(111)
    plt.plot(X, D, color="blue", linewidth=1.0, linestyle="-")

    plt.show()
#plot_attach_one()

def plot_attach():
    bta = 2.
    k = 10
    k_on  = 10.
    X = np.linspace(-3.,3., 2560)

    D = k_plus(X, 0, 2, bta, k, k_on,w=True)
    plt.figure(figsize=(8,6), dpi=80)
    plt.subplot(111)
    plt.plot(X, D, color="blue", linewidth=1.0, linestyle="-")

    plt.show()
#plot_attach()
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

#returns sth like [-n_neigh ... -1, 1 ... n_neigh+1], so the list of p taken into account
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
        raise ValueError('arg1 and arg2 must have the same length')

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
        raise ValueError('Probs should add to 1; input dicts should have the form length, values=[0., 1.], probabilities=[0.3, 0.7]')
    if not len(values) == len(probabilities):
        raise ValueError('arg1 and arg2 must have the same length')
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

def lin_fit(x, a, t):
    return a * x + t