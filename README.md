# SlidingFilamentSim

This software simulates *sliding filaments* [https://en.wikipedia.org/wiki/Sliding_filament_theory](https://en.wikipedia.org/wiki/Sliding_filament_theory)

In collaboration with [https://github.com/FWijanto/SlidingFilamentSimulation](https://github.com/FWijanto/SlidingFilamentSimulation)

## How to use:
1. **Choose operation mode** out of 
* fixed velocity drive (file run\_speed.py)
* fixed load drive (file run\_load.py)
* spring drive (spring anchor driven with fixed velocity) (file run\_spring.py)
2. **Set parameters:** this has to be done only in the above mentioned files in the _initalization block_ marked with \#
3. **Add runs** with the sets of parameters with simulation.addrun(parameters in here)
4. **Run simulations** with the commmand simulation.run(runId), plot stuff, etc
