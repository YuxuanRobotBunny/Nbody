clear;
clc;

% Canonical Sundman comparison on eccentric orbit:
%   1) fixed-step Verlet (symplectic baseline)
%   2) naive adaptive Verlet (not symplectic in general)
%   3) canonical Sundman + implicit midpoint (extended Hamiltonian)

params = struct('G', 1.0, 'm0', 1.0, 'masses', 1.0, 'dim', 2, 'eps_soft', 0.0);
params = nbody_prepare_params(params);

r0 = 1.0;
mu = params.G * params.m0;
v_circ = sqrt(mu / r0);
v_factor = 0.25;
vy0 = v_factor * v_circ;

T = orbital_period_from_initial(r0, vy0, mu);
n_periods = 20;
t_end = n_periods * T;

R0 = [r0, 0.0];
V0 = [0.0, vy0];

% (1) fixed-step Verlet
h_fix = T / 2500;
[t_fix, R_fix, V_fix, s_fix] = verlet_fixed(@nbody_acceleration, [0, t_end], R0, V0, h_fix, params);
[dE_fix, r_fix] = diagnostics_from_rv(t_fix, R_fix, V_fix, params);

% (2) naive adaptive Verlet
tol = 1e-4;
opts_ad = struct('rtol', tol, 'atol', 1e-10, 'fac_max', 2.0, 'fac_max_reject', 0.8, ...
    'h_min', 1e-8, 'h_max', 0.2);
[t_ad, R_ad, V_ad, s_ad] = verlet_adaptive_naive(@nbody_acceleration, ...
    [0, t_end], R0, V0, 0.02, tol, params, opts_ad);
[dE_ad, r_ad] = diagnostics_from_rv(t_ad, R_ad, V_ad, params);

% (3) canonical Sundman midpoint with g(r)=r
g_fun = @(R) max(norm(R(:), 2), 1e-8);
grad_g_fun = @(R) grad_norm_r(R);
h_tau = T / 1200;
opts_can = struct('newton_tol', 1e-12, 'newton_maxit', 20, 'g_min', 1e-8, 'g_max', 10);
[tau_can, t_can, R_can, V_can, s_can] = canonical_sundman_midpoint( ...
    mu, g_fun, grad_g_fun, [0, t_end], R0, V0, h_tau, opts_can);
[dE_can, r_can] = diagnostics_from_rv(t_can, R_can, V_can, params);

fprintf('Canonical Sundman test: v(0)=%.3f v_c, T=%.6f, periods=%.1f\n', ...
    v_factor, T, n_periods);
fprintf('\nSummary:\n');
print_row('Verlet fixed', max(abs(dE_fix)), min(r_fix), max(r_fix), s_fix.accel_evals);
print_row('Verlet adaptive naive', max(abs(dE_ad)), min(r_ad), max(r_ad), s_ad.accel_evals);
print_row('Canonical Sundman midpoint', max(abs(dE_can)), min(r_can), max(r_can), s_can.rhs_evals);
fprintf('Canonical constraint max |H+p_t| = %.3e\n', max(abs(s_can.constraint_hist)));
fprintf('Canonical Newton max it = %d, max residual = %.3e\n', s_can.newton_it_max, s_can.newton_res_max);

% -------------------------------------------------------------------------
% Figure 1: DeltaE
% -------------------------------------------------------------------------
figure('Name', 'Canonical Sundman Compare: DeltaE');
plot(t_fix, dE_fix, 'LineWidth', 1.4, 'DisplayName', 'Verlet fixed'); hold on;
plot(t_ad, dE_ad, 'LineWidth', 1.2, 'DisplayName', 'Verlet adaptive naive');
plot(t_can, dE_can, 'LineWidth', 1.2, 'DisplayName', 'Canonical Sundman midpoint');
grid on;
xlabel('t');
ylabel('\Delta E(t)');
title(sprintf('Energy drift over %.1f periods', n_periods));
legend('Location', 'best');

% -------------------------------------------------------------------------
% Figure 2: radius
% -------------------------------------------------------------------------
figure('Name', 'Canonical Sundman Compare: radius');
plot(t_fix, r_fix, 'LineWidth', 1.4, 'DisplayName', 'Verlet fixed'); hold on;
plot(t_ad, r_ad, 'LineWidth', 1.2, 'DisplayName', 'Verlet adaptive naive');
plot(t_can, r_can, 'LineWidth', 1.2, 'DisplayName', 'Canonical Sundman midpoint');
grid on;
xlabel('t');
ylabel('r(t)=||R(t)||_2');
title('Radius evolution');
legend('Location', 'best');

% -------------------------------------------------------------------------
% Figure 3: physical dt history
% -------------------------------------------------------------------------
figure('Name', 'Canonical Sundman Compare: dt');
plot(s_ad.t_attempt(s_ad.accepted_attempt), s_ad.h_attempt(s_ad.accepted_attempt), '.-', ...
    'DisplayName', 'Adaptive naive accepted dt'); hold on;
if ~isempty(s_can.dt_phys_hist)
    plot(t_can(1:end-1), s_can.dt_phys_hist, '.-', 'DisplayName', 'Canonical Sundman dt');
end
grid on;
xlabel('t');
ylabel('physical dt');
title('Physical time-step histories');
legend('Location', 'best');

% -------------------------------------------------------------------------
% Figure 4: canonical constraint
% -------------------------------------------------------------------------
figure('Name', 'Canonical Sundman Compare: constraint');
plot(t_can, s_can.constraint_hist, 'LineWidth', 1.3);
grid on;
xlabel('t');
ylabel('C(t) = H + p_t');
title('Extended-system constraint drift');

% Optional: t(tau) map
figure('Name', 'Canonical Sundman Compare: t vs tau');
plot(tau_can, t_can, 'LineWidth', 1.3);
grid on;
xlabel('\tau');
ylabel('t');
title('Canonical Sundman map t(\tau)');

function g = grad_norm_r(R)
R = R(:);
r = norm(R, 2);
if r <= 1e-12
    g = zeros(size(R));
else
    g = R / r;
end
end

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

function print_row(name, max_dE, rmin, rmax, evals)
fprintf('%-28s max|DeltaE|=%-10.3e  r in [%.6f, %.6f]  eval=%d\n', ...
    name, max_dE, rmin, rmax, evals);
end
