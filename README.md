
# SlidingFilamentSim

This software simulates *sliding filaments* [https://en.wikipedia.org/wiki/Sliding_filament_theory](https://en.wikipedia.org/wiki/Sliding_filament_theory)

Supervised by [https://github.com/FWijanto/SlidingFilamentSimulation](https://github.com/FWijanto/SlidingFilamentSimulation)

## How to use:

#### **Choose operation mode** in here:

 - fixed velocity run (<i class="icon-file"></i> `run_speed.py`) 
 - fixed load run (<i class="icon-file"></i> `run_load.py`) 
 - spring coupled run (spring anchor driven with fixed velocity) (<i class="icon-file"></i>
   `run_spring.py`)

> **Hint:** You can choose between different options: 
> **fixed velocity**
> - *const*: within one run of the simulation, the velocity, the filament is driven with, is constant
> - *poly*: the velocity can be defined by any polynomial with given coefficients
> **spring**
> - *stiff*: assume a spring with infinite stiffness. Therefore, the simulation converges to the velocity mode with the difference, that the calculation is exact rather than an approximation

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
```
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
```
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
```
repetitions_with_same_parameters = 10
#configure mulitprocessing TO BE USED WITH CAUTION ----> RAM OVERFLOW
n_cores = 8
use_multiprocessing = False
```

#### **Add runs** with the sets of parameters with simulation.addrun(parameters in here)
5. **Run simulations** with the commmand simulation.run(runId), plot stuff, etc

    
