# N-body 

This folder provides a clean starter implementation for the Newton/N-body part:

- Generic first-order state layout:
  - `y = [R(:); V(:)]`
  - `R, V` are `N x d` with `d=2` or `d=3`
- Physics:
  - `nbody_rhs.m`
  - `nbody_acceleration.m`
  - `nbody_energy.m`
  - `nbody_angular_momentum.m`
  - `nbody_min_distances.m`
- Integrators:
  - `bs23_step.m`
  - `bs23_fixed.m`
  - `bs23_adaptive.m`
  - `verlet_fixed.m`
  - `verlet_adaptive_naive.m`
  - `sundman_verlet_fixed.m`
  - `canonical_sundman_midpoint.m`

## Demo scripts

- `demo_single_orbit_fixed.m`
  - fixed step order test for circular orbit
  - binary search for largest step satisfying a target error
  - Richardson extrapolation comparison

- `demo_single_orbit_adaptive.m`
  - adaptive BS-RK23 for multiple eccentricities
  - tolerance sweep and accepted/rejected-step diagnostics

- `demo_adaptive_explain.m`
  - circular exact-vs-numerical comparison
  - normalized one orbit x(t) comparison
  - periapsis vs step size relation for the hardest case

- `demo_nbody_case.m`
  - `N=5` multi-body run with trajectory, energy, angular momentum diagnostics
  - optional MP4 movie output

- `demo_verlet_vs_bs23.m`
  - long time comparison of velocity-Verlet vs BS-RK23
  - trajectory, energy drift, and cost comparison

- `demo_symplectic_energy_compare.m`
  - project-A style long time energy comparison
  - fixed step symplectic vs variable step and non-symplectic methods

- `demo_phase_space_area.m`
  - project-C style phase-space area preservation experiment
  - cloud-area evolution in (q,p) for oscillator
  - includes state dependent variable step maps to show geometric drift

- `demo_sundman_verlet_compare.m`
  - high eccentricity comparison: fixed Verlet vs naive adaptive Verlet vs Sundman-Verlet
  - compares energy drift, radius drift, and effective physical step size

- `demo_canonical_sundman_compare.m`
  - extended-phase canonical Sundman midpoint prototype
  - compares fixed Verlet, naive adaptive Verlet, and canonical Sundman
  - reports constraint drift C=H+p_t

## Suggested execution order

1. `demo_single_orbit_fixed`
2. `demo_single_orbit_adaptive`
3. `demo_adaptive_explain`
4. `demo_nbody_case`
5. `demo_verlet_vs_bs23`
6. `demo_symplectic_energy_compare`
7. `demo_phase_space_area`
8. `demo_sundman_verlet_compare`
9. `demo_canonical_sundman_compare`
