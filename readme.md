# Polymer dynamics simulations
This repository contains code for performing simulations of a semi-flexible polymer:
* Monte-Carlo. It uses the Metropolis-Hasting algorithm. Configurations are generated through crankshaft and pivot moves. The elasticity is not taken into account and the bond length between monomers remains unchanged.
* Langevin. It uses Langevin dynamics. An integration time step needs to be specified.

## Compilation
Compilation of both implementations is achieved through:
```
g++ -std=c++11 main.cpp -o prog
```

## Execution
The program is then executed as follows:
```
./prog < config
```
where `config` is the configuration file. Examples of such files for both implementations are given in the `configs` folder.

This program prints the energy of the system at each time step. The configurations are written in the `pol.xyz` file.
