clear;
clc;

% Compare three Verlet variants on a high-eccentricity single-orbit problem:
%   1) fixed-step Verlet (symplectic baseline)
%   2) naive adaptive Verlet (step-doubling control)
%   3) Sundman-transformed fixed-tau Verlet
%
% Goal: show why "adaptive Verlet" may lose geometric quality, while Sundman
% can be a practical middle ground for strongly varying orbital speed.

params = struct('G', 1.0, 'm0', 1.0, 'masses', 1.0, 'dim', 2, 'eps_soft', 0.0);
params = nbody_prepare_params(params);

r0 = 1.0;
mu = params.G * params.m0;
v_circ = sqrt(mu / r0);
v_factor = 0.25;
vy0 = v_factor * v_circ;

T = orbital_period_from_initial(r0, vy0, mu);
n_periods = 40;
t_end = n_periods * T;

R0 = [r0, 0.0];
V0 = [0.0, vy0];

% Baseline fixed-step Verlet.
h_fix = T / 2500;
[t_fix, R_fix, V_fix, s_fix] = verlet_fixed(@nbody_acceleration, [0, t_end], R0, V0, h_fix, params);
[dE_fix, r_fix] = diagnostics_from_rv(t_fix, R_fix, V_fix, params);

% Naive adaptive Verlet.
tol = 1e-4;
opts_ad = struct('rtol', tol, 'atol', 1e-10, 'fac_max', 2.0, ...
    'fac_max_reject', 0.8, 'h_min', 1e-8, 'h_max', 0.2);
[t_ad, R_ad, V_ad, s_ad] = verlet_adaptive_naive(@nbody_acceleration, ...
    [0, t_end], R0, V0, 0.02, tol, params, opts_ad);
[dE_ad, r_ad] = diagnostics_from_rv(t_ad, R_ad, V_ad, params);

% Sundman fixed-tau Verlet: dt = g(r) dtau, choose g(r)=r.
h_tau = T / 300;
g_fun = @(R, ~) max(norm(R(1, :), 2), 1e-6);
opts_sd = struct('g_min', 1e-6, 'g_max', 10.0);
[tau_sd, t_sd, R_sd, V_sd, s_sd] = sundman_verlet_fixed(@nbody_acceleration, g_fun, ...
    [0, t_end], R0, V0, h_tau, params, opts_sd);
[dE_sd, r_sd] = diagnostics_from_rv(t_sd, R_sd, V_sd, params);

fprintf('High-eccentricity setup: v(0)=%.3f v_c, T=%.6f, simulated periods=%.1f\n', ...
    v_factor, T, n_periods);
fprintf('\nSummary:\n');
print_row('Verlet fixed', max(abs(dE_fix)), min(r_fix), max(r_fix), s_fix.accel_evals, h_fix, h_fix);
print_row('Verlet adaptive naive', max(abs(dE_ad)), min(r_ad), max(r_ad), s_ad.accel_evals, ...
    min(s_ad.h_accepted), max(s_ad.h_accepted));
print_row('Sundman-Verlet (g=r)', max(abs(dE_sd)), min(r_sd), max(r_sd), s_sd.accel_evals, ...
    s_sd.dt_phys_min, s_sd.dt_phys_max);
fprintf('\nInterpretation hint:\n');
fprintf('- Fixed Verlet is geometric baseline (constant dt).\n');
fprintf('- Adaptive naive Verlet improves local control but may enlarge aphelion drift.\n');
fprintf('- This Sundman prototype regularizes dt near periapsis, but is not a fully\n');
fprintf('  canonical symplectic time-transform implementation, so energy can still drift.\n');

% -------------------------------------------------------------------------
% Figure 1: energy drift
% -------------------------------------------------------------------------
figure('Name', 'Sundman vs Adaptive Verlet: DeltaE');
plot(t_fix, dE_fix, 'LineWidth', 1.4, 'DisplayName', 'Verlet fixed'); hold on;
plot(t_ad, dE_ad, 'LineWidth', 1.2, 'DisplayName', 'Verlet adaptive naive');
plot(t_sd, dE_sd, 'LineWidth', 1.2, 'DisplayName', 'Sundman-Verlet (g=r)');
grid on;
xlabel('t');
ylabel('\Delta E(t)');
title(sprintf('Energy behavior over %.1f periods', n_periods));
legend('Location', 'best');

% -------------------------------------------------------------------------
% Figure 2: radius drift
% -------------------------------------------------------------------------
figure('Name', 'Sundman vs Adaptive Verlet: radius');
plot(t_fix, r_fix, 'LineWidth', 1.4, 'DisplayName', 'Verlet fixed'); hold on;
plot(t_ad, r_ad, 'LineWidth', 1.2, 'DisplayName', 'Verlet adaptive naive');
plot(t_sd, r_sd, 'LineWidth', 1.2, 'DisplayName', 'Sundman-Verlet (g=r)');
grid on;
xlabel('t');
ylabel('r(t)=||R(t)||_2');
title('Radius drift comparison');
legend('Location', 'best');

% -------------------------------------------------------------------------
% Figure 3: physical step-size behavior
% -------------------------------------------------------------------------
figure('Name', 'Sundman vs Adaptive Verlet: physical dt');
plot(s_ad.t_attempt(s_ad.accepted_attempt), s_ad.h_attempt(s_ad.accepted_attempt), '.-', ...
    'DisplayName', 'Adaptive naive accepted dt'); hold on;
if ~isempty(s_sd.dt_phys_hist)
    plot(t_sd(1:end-1), s_sd.dt_phys_hist, '.-', 'DisplayName', 'Sundman effective dt');
end
grid on;
xlabel('t');
ylabel('physical step size dt');
title('Physical time-step history');
legend('Location', 'best');

% Optional: tau-to-t relation for Sundman.
figure('Name', 'Sundman tau-to-t map');
plot(tau_sd, t_sd, 'LineWidth', 1.4);
grid on;
xlabel('\tau');
ylabel('t');
title('Sundman map: physical time vs fictitious time');

function T = orbital_period_from_initial(r0, vy0, mu)
E = 0.5 * vy0^2 - mu / r0;
if E >= 0
    error('Initial condition is not a bound orbit (E >= 0).');
end
a = -mu / (2 * E);
T = 2 * pi * sqrt(a^3 / mu);
end

function [dE, rnorm] = diagnostics_from_rv(t, R_hist, V_hist, params)
nt = numel(t);
E = zeros(nt, 1);
rnorm = zeros(nt, 1);
for k = 1:nt
    Rk = R_hist(:, :, k);
    Vk = V_hist(:, :, k);
    E(k) = nbody_energy(pack_state(Rk, Vk), params);
    rnorm(k) = norm(Rk(1, :), 2);
end
dE = E - E(1);
end

function print_row(name, max_dE, rmin, rmax, evals, hmin, hmax)
fprintf('%-24s max|DeltaE|=%-10.3e  r in [%.6f, %.6f]  eval=%d  dt in [%.3e, %.3e]\n', ...
    name, max_dE, rmin, rmax, evals, hmin, hmax);
end


