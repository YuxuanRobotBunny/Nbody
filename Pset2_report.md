
This submission covers both parts of `Pset2.pdf`:

1. Temporal discretization of KdV with pseudo-spectral spatial discretization.
2. Newton/N-body dynamics with Bogacki-Shampine RK23 and Verlet integrators.

All computations were carried out in MATLAB. The code is organized so that the KdV and N-body solvers are reusable beyond the specific experiments required by the assignment.

## Part 1: KdV

### Spatial discretization

I discretized the periodic domain `[-L/2, L/2]` with a Fourier pseudo-spectral method. The state is stored in Fourier space as `u_hat`, with the KdV ODE written as

```text
d/dt u_hat = A u_hat + B(u_hat),
A = i k^3,
B(u_hat) = -3 i k FFT( (IFFT(u_hat))^2 ).
```

For even grid sizes, the unmatched Nyquist mode is set to zero when odd derivatives are formed. This avoids spurious evolution of that special mode.

### Time integrators

I implemented two second-order time integrators.

- `SBDF2`: an implicit-explicit multistep method using BDF2 for the linear term and AB2 for the nonlinear term.
- `ETDRK2`: an exponential RK2 method with a predictor-corrector structure.

The ETDRK2 implementation uses stable evaluations of the `phi`-functions for small `|lambda dt|`, replacing direct formulas with Taylor expansions when needed.

### Why the linear term is stiff

In Fourier space the linear term is `i k^3 u_hat`. As the number of modes increases, the largest wavenumber scales like `|k_max| ~ N/L`, so the fastest oscillation frequency grows like `|k_max|^3`. This forces explicit methods to use very small time steps to resolve high-frequency behavior even when the physical solution remains smooth. That is the stiffness associated with the dispersive linear term.

### Stability observations

I used `demo_kdv_stability.m` to empirically scan stable time steps across several resolutions. The scan is based on repeated runs over multiple soliton periods and declares instability once the numerical amplitude blows up or produces non-finite values.

The qualitative trend is that the admissible time step decreases as the number of Fourier modes increases. This is consistent with the cubic growth of the linear frequency.

### Accuracy and successive refinements

I used `demo_kdv_accuracy.m` with a single soliton and quarter-period evolution time. For each method I computed:

- successive-refinement errors `||u_dt - u_dt/2||_inf`,
- final-time errors against the exact traveling-wave solution,
- estimated orders from step halving.

The expected asymptotic behavior is second-order convergence in time for both ETDRK2 and SBDF2 when the temporal error dominates the spatial error.

### PDE accuracy and robustness

I used `demo_kdv_robustness.m` to reduce the number of modes while keeping the time step at one half of the empirical stability threshold. The final-time numerical profiles were compared directly against the exact soliton profile after one full period.

The main conclusion is that a method can look stable but still lose shape fidelity when the spatial resolution is too low. This is why robustness should be judged from full profile comparisons, not only from the absence of blow-up.

### Two-soliton interaction

I used `demo_kdv_two_solitons.m` to simulate an overtaking collision between two solitons of different speeds. The script produces snapshot figures and can optionally write an MP4 movie. The setup was chosen so that the solitons are initially well-separated and interact during the simulated time window.

## Part 2: Newton/N-body dynamics

### First-order reformulation

The equations were written as a first-order system using the state

```text
y = [R(:); V(:)].
```

The solver interface is generic in the sense that the RK23 integrator accepts an arbitrary first-order right-hand side `rhs_fun(t, y, params)`.

### Fixed-step circular orbit

I verified the order of the fixed-step BS-RK23 implementation with `demo_single_orbit_fixed.m`. The measured convergence rate is essentially third order, matching the theoretical order of the method.

I also used the exact circular orbit to determine the largest time step such that the maximum position error remains below `10^-2 r`.

### Adaptive integration

I corrected the adaptive BS-RK23 controller so that the PI step-size logic uses the previous accepted-step error instead of reusing the current one.

Then I used `demo_single_orbit_adaptive.m` to calibrate an absolute position tolerance for each eccentricity level and verify the final return-to-start error after one orbital period. This is more faithful to the assignment than directly identifying the input tolerance with the achieved global error.

### Multiple planets

The N-body experiment is implemented in `demo_nbody_case.m`. The setup uses five planets with mildly perturbed near-circular initial conditions so that the trajectories remain visually interesting without immediately causing singular close encounters.

The script reports:

- accepted and rejected adaptive steps,
- total RHS evaluations,
- minimum distance to the star,
- minimum pairwise distance,
- energy drift,
- angular-momentum drift.

It can also optionally produce an MP4 movie.

### Long-time integration and symplectic comparison

I used `demo_verlet_vs_bs23.m` to compare fixed-step Verlet and fixed-step BS-RK23 over many orbital periods for circular and elliptical cases.

The script now uses a common physical end time for all tested step sizes, making the energy-drift comparison fair. The reference step size is recomputed from the fixed-step circular-orbit target instead of being hard-coded.

The expected qualitative result is:

- Verlet shows bounded oscillatory energy error over long times,
- BS-RK23 generally has stronger secular drift at the same step size,
- despite being only second order, the symplectic method can be more reliable for long-time Hamiltonian dynamics.

## Files

### KdV

- `kdv_setup.m`
- `kdv_soliton.m`
- `kdv_nonlinear_hat.m`
- `kdv_rhs_hat.m`
- `etdrk2_kdv.m`
- `sbdf2_kdv.m`
- `demo_kdv_stability.m`
- `demo_kdv_accuracy.m`
- `demo_kdv_robustness.m`
- `demo_kdv_two_solitons.m`

### N-body

- `bs23_step.m`
- `bs23_fixed.m`
- `bs23_adaptive.m`
- `verlet_fixed.m`
- `nbody_rhs.m`
- `nbody_acceleration.m`
- `nbody_energy.m`
- `nbody_angular_momentum.m`
- `demo_single_orbit_fixed.m`
- `demo_single_orbit_adaptive.m`
- `demo_nbody_case.m`
- `demo_verlet_vs_bs23.m`

## Notes

- The KdV scripts are written to produce the figures and tables needed for the assignment, but the exact “best” time-step choices should be determined from the actual output on the target machine.
- De-aliasing is optional here and was left disabled by default so that the implementation remains close to the assignment statement.
- The optional movie-writing flags are disabled by default to keep the scripts lightweight.
