
# SlidingFilamentSimulation

This software simulates *sliding filaments* [https://en.wikipedia.org/wiki/Sliding_filament_theory](https://en.wikipedia.org/wiki/Sliding_filament_theory)

Supervised by [https://github.com/FWijanto/SlidingFilamentSimulation](https://github.com/FWijanto/SlidingFilamentSimulation)

For further information just create an issue or email me :-)

## How to use:

<<<<<<< HEAD
### **Prepare** the script

#### **Choose operation mode** by choosing the right file:

 - fixed velocity run (`run_speed.py`) 
 - fixed load run ( `run_load.py`) 
 - spring coupled run (spring anchor driven with fixed velocity) (`run_spring.py`)
=======
#### **Choose operation mode** in here:

 - fixed velocity run (<i class="icon-file"></i> `run_speed.py`) 
 - fixed load run (<i class="icon-file"></i> `run_load.py`) 
 - spring coupled run (spring anchor driven with fixed velocity) (<i class="icon-file"></i>
   `run_spring.py`)
>>>>>>> 8adc779595dd9bc12bb54b9b6015f8f553476963

> **Hint:** You can choose between different options: 
> 1. **fixed velocity**
> - *const*: within one run of the simulation, the velocity, the filament is driven with, is constant
> - *poly*: the velocity can be defined by any polynomial with given coefficients
> 2. **spring**
> - *stiff*: assume a spring with infinite stiffness. Therefore, the simulation converges to the velocity mode with the difference, that the calculation is exact rather than an approximation

<<<<<<< HEAD
####  **Rename your simulation** and specify the results path
Within the chosen document, look for this section:
```python
option = ''
name = 'please rename'
path_to_results_directory = 'please fill carefully or leave the string blank'
```
Change the option variable, if you want (see hint above). Then, choose a name for the simulation (in the results folder, a subdirectory will be created, named: date, time, name). Finally, fill in the path to the desired results directory (if you leave it blank, results will be stored in the workdir/res/)

#### **Set** the **parameters** 
Look for a section resembling the following:

```python
=======
#### <i class="icon-pencil"></i> **Rename a document**
Within the chosen document, look for this section:
```
option = ''
name = 'please rename'
path_to_results_directory = 'please fill carefully leave the string blank'
```
Change the option variable, if you want (see hint above). Then, choose a name for the simulation (in the results folder, a subdirectory will be created, named: date, time, name). Finally, fill in the path to the desired results directory (if you leave it blank, results will be stored in the workdir/res/)

#### <i class="icon-cog"></i>**Set parameters** 
Look for a section resembling the following:

```
>>>>>>> 8adc779595dd9bc12bb54b9b6015f8f553476963
#store data in ram / write them to text files?
s_store = False
p_store = True
f_store = False
sum_f_store = True
sum_p_store = False
pos_store = True
writeText = True
```
Here, setup which variables you want to store. If you just want to see plots in the end, no need to store anything on the disk: set `writeText = False`; else, all observables above set to `True` will be stored in `.dat` files.
<<<<<<< HEAD

---
```python
=======
```
>>>>>>> 8adc779595dd9bc12bb54b9b6015f8f553476963
#most important parameters
n_heads = int(1e2)
n_iterations_per_simulation = int(1e4)
delta_t = 5e-3
beta = 2.
kappa = 10.
k_on = 10.
neighborhood_criterion = 0.01
start_time = 0.
distance_between_binding_sites = 2.
random.seed(121155)
```
In here, all the physical parameters can be set. Further explanations can be found in the report.
<<<<<<< HEAD

---
```python
=======
```
>>>>>>> 8adc779595dd9bc12bb54b9b6015f8f553476963
############################################
#####CASE OPTION CONST######################
#scan for load between
step_min_val = 0
#and
step_max_val = 20
#number of steps in between
step_n_jumps = 4

#####CASE OPTION POLY#######################
coeff_of_velocity_polynomial = [[5, 0.], [10, 0.]]
############################################
```
This is the exemplary part of the `run_speed.py` file. The first part is for scanning a range of different driving velocities. The given example will return later in the code `[0, 5, 10, 15, 20]`

Alternatively, you can set the coefficients $c_i$ of the polynomial $v(t) = c_0 + c_1 \cdot t + \dots + c_n \cdot t^n$
<<<<<<< HEAD

---
```python
=======
```
>>>>>>> 8adc779595dd9bc12bb54b9b6015f8f553476963
repetitions_with_same_parameters = 10
#configure mulitprocessing TO BE USED WITH CAUTION ----> RAM OVERFLOW
n_cores = 8
use_multiprocessing = False
```

<<<<<<< HEAD
In this very last parameter part, you can define, how many repetitions you want to run (in order to get significant results!). If you feel like experimenting around a bit with multiprocessing to get the simulation done in shorter time, don't hesitate to to set it true. If you do so, make sure, you don't start to many simulations at once (memory overflows can cause your os to crash), also the `n_cores` should not exceed the number of physical cores - 1.

#### **Add runs** to the cue
There we are : )
```python
sim = simulation(   mode, ...)
if option == 'const':
        velocities = [i * (step_max_val - step_min_val) / step_n_jumps + step_min_val for i in range(step_n_jumps + 1)]
        for r in range(repetitions_with_same_parameters):
            for v in velocities:
                sim.add_run(v=v, n_sim=r)
        n = [i+1 for i in range(len(velocities) * repetitions_with_same_parameters)]
```
In this block, an object **sim** of the imported simulation class is created (of course, we have to pass all the parameters to it), then, all the runs (length of the velocity / load list times number of repetitions) will be added to **sim**.

---
```python
    if not use_multiprocessing:
        for ni in n:
            sim.start_run(ni)
    else:
        pool = Pool(processes=n_cores)
        pool.map(sim.start_run, n)
        for nr in n:
            sim.start_run(nr)
```
This will actually execute the simulation

---
```python
    ########get average forces for fixed velocity
    for i in n:
        res = sim.average_norm_force_single(i)
    sim.average_norm_force() #equivalent to above
    sim.plot_f_norm__v() #plot f_normalized(v)
    
    ########get average velocities for fixed load
    for i in n:
        average_velocity_single(i)
    sim.plot_f_norm()
    sim.plot_v__f_norm() #plot sliding_elocity(load)

    #other plotting stuff
    sim.plot_pos(n) #plot the filament position(t)
    sim.plot_p(n) #plot the number of attached heads(t)
    sim.plot_f(n) #plot the total load applied to the filament(t)
    sim.plot_f_norm #same but load is normalized by the total number of heads
    
    ########experimental part, see option "step in simulation3.py" (if you're already there and don't know, what's going on, feel free to for this in an issue)    
    sim.plot_f_v_step()
    sim.plot_v_step()
```
Last but not least: the **post processing** and **plotting**...

---
### **Run** the script
> **Important notes:**
> 1. Make sure, you have python3 including the packages `numpy, scipy, matplotlib, multiprocessing` installed 
> 2. If you want to use multiprocessing, don't run the script in an iPython console or notebook
> 3. Start with a small number of runs (memory overflow possible, you don't want to wait for your pc when it is trying to store and fetch all the state variables on the disk ;-)

Apart from the notes above, it's just opening a terminal, navigate to the workdir and type `python3 run_speed.py` or one of the other config scripts :-)
=======
#### **Add runs** with the sets of parameters with simulation.addrun(parameters in here)
5. **Run simulations** with the commmand simulation.run(runId), plot stuff, etc
>>>>>>> 8adc779595dd9bc12bb54b9b6015f8f553476963
