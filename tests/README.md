# Tests

- `make -C langevin check` builds `bin/check_forces` from `check_forces.cpp` and runs it.
  For each force field it builds a small random configuration and compares the force from
  `energy_force()` with the central difference $-[U(x+\epsilon e)-U(x-\epsilon e)]/(2\epsilon)$,
  $\epsilon = 10^{-6}$, of the energy from the same call. Error $= \max|F-F_{fd}| / \max(\max|F|, 1)$;
  pass if $< 10^{-5}$. Polar classes run in 2D; their "(frozen polarity)" lines are diagnostics
  (not counted) with the polarity vectors held fixed. It also checks that Point/Axis/Plane
  constraints hold to $10^{-12}$ after one `Lang_VVerlet` step. Exit status is nonzero on failure.
- `make -C langevin golden` runs `bin/simu_langevin` on `golden/config.yaml` in a temporary
  directory and compares `thermo.dat` and `positions_final.dat` with the stored copies
  (tolerance $5\cdot 10^{-8}$ on $|\Delta|/\max(|x_\mathrm{ref}|, 1)$, the print precision). `golden/compare.sh --update` regenerates the stored files.
