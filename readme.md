# Polymer dynamics simulations
This repository contains code for performing simulations of a semi-flexible polymer:
* Monte-Carlo. It uses the Metropolis-Hasting algorithm. Configurations are generated through crankshaft and pivot moves. The elasticity is not taken into account and the bond length between monomers remains unchanged.
* Langevin. It uses Langevin dynamics. An integration time step needs to be specified.

## Environment and compilation
The environment is managed with [pixi](https://pixi.sh): `pixi.toml` lists the compiler, `gsl`, `boost`, `yaml-cpp` and the Python packages, all from conda-forge. `pixi.lock` pins the exact versions.

Install pixi, then build the Langevin code from the repository root:
```
pixi run build
```
The first call creates the local environment in `.pixi/` and compiles `bin/simu_langevin`. The binary keeps the path to the environment libraries, so it runs without activating the environment.

Other commands:
* `pixi shell`: start a shell in the environment.
* `pixi run python analysis/analyze_trajectory.py ...`: run an analysis script.
* `pixi shell -e viz`: start a shell in the environment with PyMOL.

## Execution
The program is then executed as follows:
```
./simu_langevin [positions_init.dat] < config
```
where `config` is the configuration file, and `positions_init.dat` is an optional file containing the initial configuration.

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
