# Polymer simulations

Research code for coarse-grained simulations of bead-spring polymers:

- `langevin/`: Langevin dynamics (velocity-Verlet) of chains with bonded potentials (Gaussian, harmonic, FENE, Kratky-Porod), pair potentials (soft core, Lennard-Jones, polar pair potentials in 2D), confinement and position constraints.
- `montecarlo/`: Metropolis Monte-Carlo of a semiflexible chain with crankshaft and pivot moves.
- `analysis/`, `pymol/`: Python scripts to analyze and visualize trajectories.

Units: $k_B = 1$; the friction $\gamma$ is a rate (friction force $-m\gamma v$).

## Build

The environment is managed with [pixi](https://pixi.sh). From the repository root:
```
pixi run build       # bin/simu_langevin
pixi run build-mc    # bin/simu_montecarlo
pixi run check       # force/energy consistency tests
```

## Run

```
bin/simu_langevin [initial_positions.dat|.xyz] < configs/config_langevin.yaml
bin/simu_montecarlo < configs/config_montecarlo
```
The configuration files document their keys. The Langevin code writes, in the current directory:
- `thermo.dat`: iteration, number of particles, potential and kinetic energy;
- `trajectory/xyz/`, `trajectory/dat/`: positions (the directory is cleared at start);
- `trajectory/polarity/`, `trajectory/neighbor/`: polarity vectors and neighbor lists;
- `positions_final.dat`, `thermo_final.dat`.

The Monte-Carlo code writes `pol.xyz`.

Analysis, for example a contact map: `pixi run python analysis/analyze_trajectory.py trajectory/xyz --compute_contacts`.

## Notes

- Polar pair potentials (`PolarPairLJ`, `PolarPair48`, `PolarPair48_2site`) are 2D only. The polarity of a bead is the in-plane normal to its chain, $n = \hat z \times \hat t$. Before November 2026, their forces omitted the dependence of $n$ on the positions of neighboring beads, and the FENE bond force had a wrong $1/r$ factor; results obtained with earlier versions are affected.
- Cutoffs: `rc` is in units of $\sigma$ for `PairLJ` and `PolarPairLJ`, and absolute for the `*48` potentials.
- The thermal noise is uniform with the variance required by the fluctuation-dissipation relation. Constraints are applied as projections of positions and velocities.

## License

MIT, see `LICENSE`.
