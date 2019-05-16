# Polymer dynamics simulations
This repository contains code for performing simulations of a semi-flexible polymer:
* Monte-Carlo. It uses the Metropolis-Hasting algorithm. Configurations are generated through crankshaft and pivot moves. The elasticity is not taken into account and the bond length between monomers remains unchanged.
* Langevin. It uses Langevin dynamics. An integration time step needs to be specified.

## Compilation
Compilation of both implementations is achieved through:
```
make -C langevin/
```

## Execution
The program is then executed as follows:
```
./prog < config
```
where `config` is the configuration file. Examples of such files for both implementations are given in the `configs` folder.

This program prints the energy of the system at each time step. The configurations are written in the `pol.xyz` file.

# Visualization of trajectories
## Pymol
With configurations stored in different files, the following commands need to be run in the Pymol interpreter:
```
reinitialize
import glob
lst = glob.glob("path/to/xyz/state*.xyz")
lst.sort()
for s in lst: cmd.load(s)
join_states trajectory, state*, 0
delete state*
hide everything, all
alter all, vdw=0.5
rebuild
show spheres, all
```
